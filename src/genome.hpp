/*
 * Copyright [2018] [XiaoFei Zhao]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SIGNATURE_HPP
#define SIGNATURE_HPP

// #include "libpopcnt-2.2/libpopcnt.h" // the libpopcnt code that is used in this file is copied into this file

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define BITATPOS(x, pos) ((x & (1ULL << pos)) >> pos)
#define NBITS(x) (8*sizeof(x))
#define ROUNDDIV(a, b) (((a) + (b)/2) / (b))


// Start of macros and method copied from https://github.com/kimwalisch/libpopcnt

#ifdef __GNUC__
	#define GNUC_PREREQ(x, y) \
			(__GNUC__ > x || (__GNUC__ == x && __GNUC_MINOR__ >= y))
#else
	#define GNUC_PREREQ(x, y) 0
#endif

#ifndef __has_builtin
	#define __has_builtin(x) 0
#endif

/*
 * This uses fewer arithmetic operations than any other known
 * implementation on machines with fast multiplication.
 * It uses 12 arithmetic operations, one of which is a multiply.
 * http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
 */
static inline uint64_t popcount64(uint64_t x)
{
	uint64_t m1 = 0x5555555555555555ll;
	uint64_t m2 = 0x3333333333333333ll;
	uint64_t m4 = 0x0F0F0F0F0F0F0F0Fll;
	uint64_t h01 = 0x0101010101010101ll;

	x -= (x >> 1) & m1;
	x = (x & m2) + ((x >> 2) & m2);
	x = (x + (x >> 4)) & m4;

	return (x * h01) >> 56;
}

// End of macros and method copied from https://github.com/kimwalisch/libpopcnt

class Entity {
public:
	std::string name;
	std::vector<uint64_t> usigs;
	double matchprob;
	Entity() {} 
	Entity(std::string infname);
	void write(std::ofstream &binfile, std::ofstream &txtfile);
	void read(std::ifstream &binfile, std::ifstream &txtfile);
};

void print_entity(const Entity &entity) {
	std::cerr << entity.name;
	for (auto u : entity.usigs) {
		std::cerr << "\t" << u;
	}
	
	std::cerr << "\t" << entity.matchprob << "\n";
}

bool operator< (const Entity &a, const Entity &b) { return a.matchprob < b.matchprob; }

void Entity::write(std::ofstream &binfile, std::ofstream &txtfile) {
	// std::cerr << "filename = " << name << std::endl; int i = 1;
	for (auto usig : usigs) {
		// std::cerr << "\t" << std::hex << usig; if (i % 4 == 0) {std::cerr << std::endl;}; i++;
		binfile.write((char*)&usig, sizeof(uint64_t));
	}
	double matchprob2 = matchprob;
	if (matchprob2 < DBL_EPSILON) { matchprob2 = DBL_EPSILON; }
	if (matchprob2 > 1.0 - DBL_EPSILON) { matchprob2 = 1.0 - DBL_EPSILON; }

	txtfile << name << "\t" << std::scientific << matchprob2 << "\t";
	txtfile << usigs.size() << "\n";
}

void save_entities(const std::vector<Entity> &entities, const std::string &outpname) {
	std::ofstream binfile(outpname + ".dat", std::ios::out | std::ios::binary);
	std::ofstream txtfile(outpname + ".txt", std::ios::out);
	if (!binfile) {
		std::cerr << "Error: cannot open the file " << outpname + ".dat" << " for writing\n";
		exit(-1);
	}
	if (!txtfile) {
		std::cerr << "Error: cannot open the file " << outpname + ".txt" << " for writing\n";
		exit(-2);
	}
	txtfile << entities.size() << "\n";
	for (auto entity : entities) {
		entity.write(binfile, txtfile);
	}
	binfile.close();
	txtfile.close();
}

void Entity::read(std::ifstream &binfile, std::ifstream &txtfile) {
	size_t vsize;
	txtfile >> name >> matchprob >> vsize;
	//std::cout << name << "\t" << matchprob << "\t" << vsize << std::endl;
	usigs = std::vector<uint64_t>(vsize);	
	for (size_t i = 0; i < vsize; i++) {
		binfile.read((char*)&usigs[i], sizeof(uint64_t));
		//std::cout << "\tRead " << usigs[i] << std::endl;
	}
}

void load_entities(std::vector<Entity> &entities, const std::string &inpname) {
	std::ifstream binfile(inpname + ".dat", std::ios::in | std::ios::binary);
	std::ifstream txtfile(inpname + ".txt", std::ios::in);
	if (!binfile) {
		std::cerr << "Error: cannot open the file " << inpname + ".dat" << " for reading\n";
		exit(-1);
	}
	if (!txtfile) {
		std::cerr << "Error: cannot open the file " << inpname + ".txt" << " for reading\n";
		exit(-2);
	}
	size_t vsize;
	txtfile >> vsize;
	entities.clear();
	entities.reserve(vsize);
	for (size_t i = 0; i < vsize; i++) {
		Entity entity;		
		entity.read(binfile, txtfile);
		entities.push_back(entity);
	}
	binfile.close();
	txtfile.close();
}

void fillusigs(Entity &entity, const std::vector<uint64_t> &signs, size_t bbits) {
	for (size_t signidx = 0; signidx < signs.size(); signidx++) {
		uint64_t sign = signs[signidx];
		int leftshift = (signidx % NBITS(uint64_t));
		for (size_t i = 0; i < bbits; i++) {
			uint64_t orval = (BITATPOS(sign, i) << leftshift);
			entity.usigs[signidx/NBITS(uint64_t) * bbits + i] |= orval;
		}
	}
}

size_t calc_intersize0(size_t &intersize, size_t &unionsize, 
		const Entity &e1, const Entity &e2, size_t sketchsize64) {
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int samebits = 0;
	unsigned int diffbits = 0;
	while (i < e1.usigs.size() && j < e2.usigs.size()) {
		if (e1.usigs[i] == e2.usigs[j]) {
			samebits++;
			i++;j++;
		} else {
			diffbits++;
			if (e1.usigs[i] < e2.usigs[j]) {
				j++;
			} else {
				i++;
			}
		}
	}
	if (i < e1.usigs.size()) {
		diffbits += e1.usigs.size() - i;
	} else if (j < e2.usigs.size()) {
		diffbits += e2.usigs.size() - j;
	}
	intersize = samebits;
	unionsize = samebits + diffbits;
	assert(intersize + unionsize == e1.usigs.size() + e2.usigs.size());
	return ROUNDDIV((sketchsize64 * NBITS(uint64_t)) * samebits, samebits + diffbits);
}

template <class T>
T non_neg_minus(T a, T b) {
	return a > b ? (a - b) : 0;
}

const size_t calc_intersize12(const Entity &e1, const Entity &e2, const size_t sketchsize64, const size_t bbits) {
	// assert (e1.usigs.size() == e2.usigs.size());	
	// assert (e1.usigs.size() == sketchsize64 * bbits);
	size_t samebits = 0;
	for (size_t i = 0; i < sketchsize64; i++) {
		uint64_t bits = ~((uint64_t)0ULL);
		// std::cout << "bits = " << std::hex << bits << std::endl;
		for (size_t j = 0; j < bbits; j++) {
			// assert(e1.usigs.size() > i * bbits + j || !fprintf(stderr, "i=%lu j=%lu bbits=%lu vsize=%lu\n", i, j, bbits, e1.usigs.size()));
			bits &= ~(e1.usigs[i * bbits + j] ^ e2.usigs[i * bbits + j]);
			// std::cout << " bits = " << std::hex << bits << std::endl;
		}

		// The popcnt method from "libpopcnt-2.2/libpopcnt.h" gives compiling error on MacOS. 
		// samebits += popcnt(&bits, sizeof(uint64_t)); // uses  __builtin_popcountll(bits) on most architectures;

#if GNUC_PREREQ(4, 2) || __has_builtin(__builtin_popcountll) 
		samebits += __builtin_popcountll(bits);
#else
		samebits += popcount64(bits)
#endif
	}
	// std::cout << " samebits = " << std::hex << samebits << std::endl;
	const size_t maxnbits = sketchsize64 * NBITS(uint64_t); 
	const size_t expected_samebits = (maxnbits >> bbits);
	if (expected_samebits) {
		return samebits;
	}
	size_t ret = non_neg_minus(samebits, expected_samebits);
	return ret * maxnbits / (maxnbits - expected_samebits);
	
#if 0 // This commented-out code works but has lots of complexity that is not justified so far
	// return sketchsize64 * NBITS(uint64_t) - samebits;
	size_t diffbits = sketchsize64 * NBITS(uint64_t) - samebits;
	size_t multiplier = (1ULL << bbits);
	size_t divider = multiplier - 1;
	diffbits = ROUNDDIV(diffbits * multiplier, divider);
	return non_neg_minus((sketchsize64 * NBITS(uint64_t)), diffbits);
#endif
}

const double intersize_to_jaccard(const unsigned int intersize, size_t sketchsize64) {
	return (double)intersize / (double)(NBITS(uint64_t) * sketchsize64);
}

const double intersize_to_dice(const unsigned int intersize, size_t sketchsize64, double smoothfactor = (DBL_MIN*(1ULL<<30ULL))) {
	return (double)(2 * intersize + smoothfactor) / (double)(NBITS(uint64_t) * sketchsize64 + intersize + smoothfactor);
}

void intersize_to_mutdist_init(std::vector<double> &intersize_to_mutdist, size_t sketchsize64, size_t kmerlen) {
	intersize_to_mutdist = std::vector<double>(NBITS(uint64_t) * sketchsize64 + 1);
	for (unsigned int i = 0; i < NBITS(uint64_t) * sketchsize64 + 1; i++) {
		double dice = intersize_to_dice(i, sketchsize64);
		double mutdist = -1.0 / kmerlen * log(dice);
		intersize_to_mutdist[i] = mutdist;
	}
	intersize_to_mutdist[NBITS(uint64_t) * sketchsize64] = 0;
}

#endif
