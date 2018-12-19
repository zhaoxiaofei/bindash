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

#ifndef HASHUTILS_HPP
#define HASHUTILS_HPP

#include "rollinghashcpp/cyclichash.h"

#include "cbuffer.hpp"

#include <algorithm>
#include <queue>
#include <set>
#include <vector>

#define MEDIAN2X(data, n) ((data[(n-1)/2] + data[n/2]))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define NBITS(x) (8*sizeof(x))

//#define ORDER(a, b) if ((a) > (b)) { (a) = (a)^(b); (b) = (a)^(b); (a) = (a)^(b); }

const uint64_t SIGN_MOD = (1ULL << 61ULL) - 1ULL; //1000000007; // (1ULL << 48ULL); // (1ULL << 49ULL) - 1ULL; // << (1L << 31L) - 1L; 

inline uint64_t doublehash(uint64_t hash1, uint64_t hash2) { return (hash1 + hash2) % SIGN_MOD; }

void binsign(std::vector<uint64_t> &signs, const uint64_t sign, const uint64_t binsize) {
	// std::cerr << "binning " << sign << std::endl;
	uint64_t binidx = sign / binsize;
	// assert(binidx < signs.size());
#ifdef CANONICAL_DENSIFICATION
	while (sign < signs[binidx]) {
		signs[binidx] = sign;
		if (0 == binidx) { break; }
		binidx--;
	}
#else
	signs[binidx] = MIN(signs[binidx], sign);
#endif
}

uint64_t univhash2(uint64_t s, uint64_t t) {
	uint64_t x = (1009) * s + (1000*1000+3) * t;
	return (48271 * x + 11) % ((1ULL << 31) - 1);
}


/* This densitfication strategy is found at
@article{shrivastava2017optimal,
  title={Optimal densification for fast and accurate minwise hashing},
  author={Shrivastava, Anshumali},
  journal={arXiv preprint arXiv:1703.04664},
  year={2017}
}
*/
int densifybin(std::vector<uint64_t> &signs) {
	uint64_t minval = UINT64_MAX;
	uint64_t maxval = 0;
	for (auto sign : signs) { 
		minval = MIN(minval, sign);
		maxval = MAX(maxval, sign);
	}
	if (UINT64_MAX != maxval) { return 0; }
	if (UINT64_MAX == minval) { return -1; }
	for (uint64_t i = 0; i < signs.size(); i++) {
		uint64_t j = i;
		uint64_t nattempts = 0;
		while (UINT64_MAX == signs[j]) {
			j = univhash2(i, nattempts) % signs.size();
			nattempts++;
		}
		signs[i] = signs[j];
	}
	return 1;
}

void ppush(std::priority_queue<uint64_t> &signqueue, std::set<uint64_t> &signset, uint64_t sign, const size_t sketchsize64) {
	if (signqueue.size() < sketchsize64 * NBITS(uint64_t) && signset.find(sign) == signset.end()) {
		signset.insert(sign);
		signqueue.push(sign);
	} else if (signqueue.top() > sign && signset.find(sign) == signset.end()) {
		signset.erase(signqueue.top());
		signqueue.pop();
		signset.insert(sign);
		signqueue.push(sign);
	}
}

const uint64_t ipow(uint64_t base, uint64_t exp) {
	uint64_t result = 1;
	while (exp) {
		if (exp & 1) { result *= base; }
		exp >>= 1;
		base *= base;
	}
	return result;
}

class _DnaPerfectHash {
public:
	uint64_t hashvalue = 0;
	const size_t kmerlen;	
	const uint64_t maxnomval;
	uint64_t NUC_TO_ID[256] = {0};

	_DnaPerfectHash(size_t k) : kmerlen(k), maxnomval(ipow(5, k-1)) {
		NUC_TO_ID['A'] = NUC_TO_ID['a'] = 1;
		NUC_TO_ID['C'] = NUC_TO_ID['c'] = 2;
		NUC_TO_ID['G'] = NUC_TO_ID['g'] = 3;
		NUC_TO_ID['T'] = NUC_TO_ID['t'] = 4;
		// NUC_TO_ID['U'] = NUC_TO_ID['u'] = 4; // RNA exluded
	}

	void eat(unsigned char);
	void update(unsigned char, unsigned char);
	void reverse_update(unsigned char, unsigned char);
};

void _DnaPerfectHash::eat(unsigned char nucleotide) {
	hashvalue *= 5;
	auto nucid = NUC_TO_ID[nucleotide];
	assert(nucid < 5);
	hashvalue += nucid;
}

void _DnaPerfectHash::update(unsigned char oldnuc, unsigned char newnuc) {
	uint64_t oldid = NUC_TO_ID[oldnuc];
	uint64_t newid = NUC_TO_ID[newnuc];
	hashvalue -= oldid * maxnomval;
	hashvalue *= 5;
	hashvalue += newid;
}

void _DnaPerfectHash::reverse_update(unsigned char newnuc, unsigned char oldnuc) {
	uint64_t newid = NUC_TO_ID[newnuc];
	hashvalue /= 5;
	hashvalue += newid * maxnomval;
}

void hashinitDna(CBuf & cbuf, _DnaPerfectHash &hf, _DnaPerfectHash &hfrc, size_t kmerlen) {
	hf.hashvalue = 0;
	hfrc.hashvalue = 0;
	for (size_t k = 0; k < kmerlen; ++k) {
		hf.eat(cbuf.getith(k));
		hfrc.eat(RCMAP[(int)cbuf.getith(kmerlen - k - 1)]);
	}
}

template <class T>
void hashupdateDna(CBuf &cbuf, std::set<uint64_t, T> &signs, 
		_DnaPerfectHash &hf, _DnaPerfectHash &hfrc, bool isstrandpreserved) {
	hf.update(cbuf.getout(), cbuf.getnewest());
	hfrc.reverse_update(RCMAP[(int)cbuf.getnewest()], RCMAP[(int)cbuf.getout()]);
	if (cbuf.slen >= cbuf.size) {
		auto signval = hf.hashvalue;
		if (!isstrandpreserved) {
			auto signval2 = hfrc.hashvalue;
			//std::cerr << "kmer = " << cbuf.tostring() << " value = " << signval << "\t" << signval2 << std::endl;
			signval = MIN(signval, signval2); // uniformity is not needed because this is a perfect hash function
		}
		signs.insert(signval);
	}
}


void hashinit0(CBuf &cbuf, CyclicHash<uint64_t> &hf, CyclicHash<uint64_t> &hfrc, size_t kmerlen) {
	hf.hashvalue = 0;
	hfrc.hashvalue = 0;
	for (size_t k = 0; k < kmerlen; ++k) {
		hf.eat(cbuf.getith(k));
		hfrc.eat(RCMAP[(int)cbuf.getith(kmerlen - k - 1)]);
	}
}

void hashupdate0(CBuf &cbuf, std::priority_queue<uint64_t> &signqueue, std::set<uint64_t> &signset, 
		CyclicHash<uint64_t> &hf, CyclicHash<uint64_t> &hfrc, 
		bool isstrandpreserved, size_t sketchsize64) {
	hf.update(cbuf.getout(), cbuf.getnewest());
	hfrc.reverse_update(RCMAP[(int)cbuf.getnewest()], RCMAP[(int)cbuf.getout()]);
	if (cbuf.slen >= cbuf.size) {
		auto signval = hf.hashvalue % SIGN_MOD;
		if (!isstrandpreserved) {
			auto signval2 = hfrc.hashvalue % SIGN_MOD;
			signval = doublehash(signval, signval2);
		}
		ppush(signqueue, signset, signval, sketchsize64);
		// std::cerr << "kmer = " << cbuf.tostring() << " value = " << signval << std::endl;
	}
}

const uint64_t estimate_genome_size0(uint64_t maxsign, size_t nsigns) {
	return SIGN_MOD / maxsign * nsigns;
}

void hashinit1(CBuf & cbuf, std::vector<CyclicHash<uint64_t>> &hfs, std::vector<CyclicHash<uint64_t>> &hfrcs, 
		size_t kmerlen, size_t sketchsize64) {
	for (size_t inc = 0; inc < sketchsize64 * NBITS(uint64_t); inc++) {	
		hfs[inc].hashvalue = 0;
		hfrcs[inc].hashvalue = 0;
		for (size_t k = 0; k < kmerlen; ++k) {
			hfs[inc].eat(cbuf.getith(k));
			hfrcs[inc].eat(cbuf.getith(kmerlen - k - 1));
		}
	}
}

void hashupdate1(CBuf & cbuf, std::vector<uint64_t> &signs, 
		std::vector<CyclicHash<uint64_t>> &hfs, std::vector<CyclicHash<uint64_t>> &hfrcs, 
		bool isstrandpreserved, size_t sketchsize64) {
	for (size_t inc = 0; inc < sketchsize64 * NBITS(uint64_t); inc++) {
		hfs[inc].update(cbuf.getout(), cbuf.getnewest());
		hfrcs[inc].reverse_update(RCMAP[(int)cbuf.getnewest()], RCMAP[(int)cbuf.getout()]);
		if (cbuf.slen >= cbuf.size) {
			auto signval = hfs[inc].hashvalue % SIGN_MOD;
			if (!isstrandpreserved) {
				auto signval2 = hfrcs[inc].hashvalue % SIGN_MOD;
				signval = doublehash(signval, signval2);
			}
			signs[inc] = MIN(signs[inc], signval);
		}
	}
}

const uint64_t estimate_genome_size1(const std::vector<uint64_t> &signs) {
	std::vector<uint64_t> deltas(signs);
	std::sort(deltas.begin(), deltas.end());
	return SIGN_MOD / MEDIAN2X(deltas, deltas.size());
}

void hashinit2(CBuf & cbuf, CyclicHash<uint64_t> &hf, CyclicHash<uint64_t> &hfrc, size_t kmerlen) {
	hf.hashvalue = 0;
	hfrc.hashvalue = 0;
	for (size_t k = 0; k < kmerlen; ++k) {
		hf.eat(cbuf.getith(k));
		hfrc.eat(RCMAP[(int)cbuf.getith(kmerlen - k - 1)]);
	}
}

void hashupdate2(CBuf & cbuf, std::vector<uint64_t> &signs, 
		CyclicHash<uint64_t> &hf, CyclicHash<uint64_t> &hfrc, 
		bool isstranspreserved, uint64_t binsize) {
	hf.update(cbuf.getout(), cbuf.getnewest());
	hfrc.reverse_update(RCMAP[(int)cbuf.getnewest()], RCMAP[(int)cbuf.getout()]);
	if (cbuf.slen >= cbuf.size) {
		// std::cerr << cbuf.tostring() << std::endl;
		auto signval = hf.hashvalue % SIGN_MOD;
		if (!isstranspreserved) {
			auto signval2 = hfrc.hashvalue % SIGN_MOD;
			signval = doublehash(signval, signval2);
		}
		binsign(signs, signval, binsize);
	}
}

const uint64_t estimate_genome_size2(const std::vector<uint64_t> &signs, size_t sketchsize64) {
	uint64_t nbins = (sketchsize64 * NBITS(uint64_t));
	uint64_t binsize = (SIGN_MOD + nbins - 1ULL) / nbins;
	std::vector<uint64_t> deltas;
	deltas.reserve(signs.size());
	uint64_t n_empty_bins = 0;
	for (uint64_t binmin = 0, i = 0; 
		binmin < SIGN_MOD && i < signs.size(); 
		binmin += binsize, i++) {
		// this is to prevent the issue https://github.com/zhaoxiaofei/bindash/issues/3 from happening
		if (signs[i] >= binmin && signs[i] < binmin + binsize) {
			deltas.push_back(signs[i] - binmin);
		} else {
			// this should rarely happens
			n_empty_bins += 1;
		}
	}
	// this is to address https://github.com/zhaoxiaofei/bindash/issues/3
	// uses linear extrapolation to compute the difference between the two consecutive hash values that span the otherwise empty bin.
	for (uint64_t i = 0; i < n_empty_bins; i++) {
		deltas.push_back(binsize * (1 + i) / (1 + signs.size() - n_empty_bins) + binsize);
	}
	assert (deltas.size() == signs.size() /*|| !fprintf(stderr, "%u != %u", deltas.size(), signs.size())*/ );
	std::sort(deltas.begin(), deltas.end());
	return SIGN_MOD / MAX(1, MEDIAN2X(deltas, deltas.size()));
}

#endif
