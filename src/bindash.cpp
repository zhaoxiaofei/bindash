/*
 * Copyright [2023] [XiaoFei Zhao] and [Jianshu Zhao]
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

#include "bhmath.hpp"
#include "cbuffer.hpp"
#include "commands.hpp"
#include "genome.hpp"
#include "hashutils.hpp"
#include "kmerset.hpp"

#include "rollinghashcpp/cyclichash.h"

//#include "kseq.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>
#include <array>
#include <queue>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define VERSION "2.2.0"

/* TODO:
 * Potential improvement to be made: CACHE_AWARE nearest neighbor search (may involve substantial additional coding).
 */
#define CACHE_SIZE (1024) // approximated
#define DIST_IS_TEMPLATED 0
#define ISIZE  (2 * CACHE_SIZE)
#define LINELEN_MAX (1024)

void *xmalloc(size_t n) {
	void *ret = malloc(n);
	if (NULL == ret) {
		fprintf(stderr, "Malloc %lu bytes failed!", n);
		abort();
	}
	return ret;
}

void *xrealloc(void *data, size_t n) {
	void *ret = realloc(data, n);
	if (NULL == ret) {
		fprintf(stderr, "Realloc %lu bytes from address %p failed!", n, data);
		abort();
	}
	return ret;
}

template <class T>
bool ispowof2(T x) { return 0 == (x & (x-1)); }

void entity_init(Entity &entity, CBuf &cbuf, const std::string &entityname, const Sketchargs &args, const size_t nseqs = SIZE_MAX) {
	entity.name = entityname;
	uint64_t genome_size = 0;
	std::minstd_rand0 g1(args.randseed);
	auto s1 = g1();
	auto s2 = g1();

	if (-1 == args.minhashtype) { // special perfect hash function for nucleotides
		std::set<uint64_t, std::greater<uint64_t>> signset;
		_DnaPerfectHash hf(args.kmerlen);
		_DnaPerfectHash hfrc(args.kmerlen);
		while (cbuf.nseqs() < nseqs) {
			cbuf.eatnext();
			if (cbuf.ceof()) { break; }
			hashupdateDna(cbuf, signset, hf, hfrc, args.isstrandpreserved);
		}
		genome_size = signset.size();
		entity.usigs = std::vector<uint64_t>(signset.begin(), signset.end());
	} else if (0 == args.minhashtype) { // multi minvals one hash function
		std::priority_queue<uint64_t> signqueue;
		std::set<uint64_t> signset;
		CyclicHash<uint64_t> hf(args.kmerlen, s1, s2, 64);
		CyclicHash<uint64_t> hfrc(args.kmerlen, s1, s2, 64);
		hashinit0(cbuf, hf, hfrc, args.kmerlen);
		while (cbuf.nseqs() < nseqs)  {
			cbuf.eatnext();
			if (cbuf.ceof()) { break; }
			hashupdate0(cbuf, signqueue, signset, hf, hfrc, args.isstrandpreserved, args.sketchsize64);
		}
		assert (signqueue.size() <= args.sketchsize64 * NBITS(uint64_t));
		if (signqueue.size() < args.sketchsize64 * NBITS(uint64_t)) {
			genome_size = signqueue.size();
		} else if (0 < signqueue.size()) {
			genome_size = estimate_genome_size0(signqueue.top(), signqueue.size());
		}
		size_t inisize = signqueue.size();
		entity.usigs = std::vector<uint64_t>(inisize, 0);
		for (size_t i = 0; i < inisize; i++) {
			entity.usigs[i] = signqueue.top();
			// std::cerr << "File" << entity.name << " stored value " << entity.usigs[i] << std::endl;
			signqueue.pop();
		}
	} else if (1 == args.minhashtype) { // one minval multi hash functions
		entity.usigs = std::vector<uint64_t>(args.sketchsize64 * args.bbits, 0);
		std::vector<uint64_t> signs(args.sketchsize64 * NBITS(uint64_t)); // carry over
		std::fill(signs.begin(), signs.end(), UINT64_MAX);
		std::vector<CyclicHash<uint64_t>> hfs;
		std::vector<CyclicHash<uint64_t>> hfrcs;
		for (size_t i = 0; i < args.sketchsize64 * NBITS(uint64_t); i++) {
			auto s1 = g1();
			auto s2 = g1();
			CyclicHash<uint64_t> hf(args.kmerlen, s1, s2, 64);
			CyclicHash<uint64_t> hfrc(args.kmerlen, s1, s2, 64);
			hfs.push_back(hf);
			hfrcs.push_back(hfrc);
		}
		hashinit1(cbuf, hfs, hfrcs, args.kmerlen, args.sketchsize64);
		while (cbuf.nseqs() < nseqs) {
			cbuf.eatnext();
			if (cbuf.ceof()) { break; }
			hashupdate1(cbuf, signs, hfs, hfrcs, args.isstrandpreserved, args.sketchsize64);
		}
		genome_size = estimate_genome_size1(signs);
		fillusigs(entity, signs, args.bbits);
	} else if (2 == args.minhashtype) { // multi binvals one hash function
		const uint64_t nbins = args.sketchsize64 * NBITS(uint64_t);
		const uint64_t binsize = (SIGN_MOD + nbins - 1ULL) / nbins;
		entity.usigs = std::vector<uint64_t>(args.sketchsize64 * args.bbits, 0);
		std::vector<uint64_t> signs(args.sketchsize64 * NBITS(uint64_t), UINT64_MAX); // carry over
		CyclicHash<uint64_t> hf(args.kmerlen, s1, s2, 64);
		CyclicHash<uint64_t> hfrc(args.kmerlen, s1, s2, 64);
		hashinit2(cbuf, hf, hfrc, args.kmerlen);
		while (cbuf.nseqs() < nseqs) {
			cbuf.eatnext();
			if (cbuf.ceof()) { break; }
			hashupdate2(cbuf, signs, hf, hfrc, args.isstrandpreserved, binsize);
		}
		// Declaration of res as an integer
		int res = 0;
		if (1 == args.dens){
			// use optimal densification
			res = opt_densify(signs);
			if (res != 0) {
				std::cerr << "Warning: the genome " << entityname << " is densified with flag " << args.dens <<  std::endl;
			}
		} else if (2 == args.dens){
			// use reverse optimal densification, or faster densification
			res = revopt_densify(signs);
			if (res != 0) {
				std::cerr << "Warning: the genome " << entityname << " is densified with flag " << args.dens <<  std::endl;
			}
		} else if (3 == args.dens){
			// use optimal densification with re-randomization
			res = rerand_densify(signs);
			if (res != 0) {
				std::cerr << "Warning: the genome " << entityname << " is densified with flag " << args.dens <<  std::endl;
			}
		}

		genome_size = estimate_genome_size2(signs, args.sketchsize64);
		fillusigs(entity, signs, args.bbits);
	}
	double entropy = bhmath_calc_entropy(cbuf.chfreqs, 256);
	entity.matchprob = bhmath_matchprob(args.kmerlen, entropy, genome_size + 1);
	// fprintf(stderr, "Entropy = %f, genome_size = %lu, matchprob = %f\n", entropy, genome_size + 1, entity.matchprob);
}

void entities_init(std::vector<Entity> &entities, std::string fname, const Sketchargs &args,
		const std::vector<std::pair<size_t, size_t>> &entityid_to_count_vec,
		const std::vector<std::string> &entityid_to_name) {
	// std::cerr << "Initiating fname " << fname << std::endl;
	CBuf cbuf(fname, args.kmerlen, args.iscasepreserved);
	size_t tot_nseqs = 0;
	for (auto entityid_to_count :  entityid_to_count_vec) {
		size_t entityid = entityid_to_count.first;
		size_t count = entityid_to_count.second;
		assert(count > 0);
		tot_nseqs += count;
		entity_init(entities[entityid], cbuf, entityid_to_name[entityid], args, tot_nseqs + 1);
	}
	if (tot_nseqs != cbuf.nseqs() && tot_nseqs != SIZE_MAX - 1) {
		std::cerr << "Runtime error: file " << fname << " expects " << tot_nseqs << " sequences but actually has " << cbuf.nseqs() << " sequences!" << std::endl;
		exit(-128);
	}
}

int cmddist_filter(double &mutdist, double &pvalue,
		size_t &raw_intersize, size_t &raw_unionsize,
		const Entity &query, const Entity &target,
		const Sketchargs &args1, const Distargs &args2,
		const std::vector<double> &intersize_to_mutdist) {

	size_t intersize;
	if (args1.minhashtype <= 0) {
		intersize = calc_intersize0(raw_intersize, raw_unionsize, query, target, args1.sketchsize64);
	} else {
		intersize = calc_intersize12(query, target, args1.sketchsize64, args1.bbits);
		raw_unionsize = NBITS(uint64_t) * args1.sketchsize64;
		if (intersize < args2.ithres) {
			raw_unionsize -= intersize;
			intersize = 0;
		}
		raw_intersize = intersize;
	}

	mutdist = intersize_to_mutdist[intersize];
	if (mutdist > args2.mthres) { return 1; }
	if (0 < intersize) {
		double andprob = query.matchprob * target.matchprob;
		double p = andprob / (query.matchprob + target.matchprob - andprob);
		pvalue = bhmath_pbinom_tail(raw_intersize, raw_unionsize, p);
	} else {
		pvalue = 1.0;
	}
	if (pvalue > args2.pthres) { return 2; }

	return 0;
}

size_t cmddist_buffer(char *&buffer, size_t &buffer_size, size_t &buffer_maxsize,
		const Entity &query, const Entity &target,
		double mutdist, const double pvalue, size_t raw_intersize, size_t raw_unionsize) {
	assert (buffer_maxsize >= LINELEN_MAX);
	if (buffer_size + LINELEN_MAX >= buffer_maxsize) {
		buffer_maxsize *= 2;
		buffer = (char*)xrealloc(buffer, buffer_maxsize);
	}
	int inc = snprintf(&buffer[buffer_size], LINELEN_MAX, "%s\t%s\t%.4e\t%.4e\t%lu/%lu\n",
			query.name.c_str(), target.name.c_str(), mutdist, pvalue, raw_intersize, raw_unionsize);
	assert(inc < LINELEN_MAX);
	assert(inc > 0);
	buffer_size += inc;
	return inc;
}

void cmddist_write(FILE *outfile, char **buffers, size_t *buffer_sizes, size_t nbuffers) {
	for (size_t i = 0; i < nbuffers; i++) {
		size_t nchars = fwrite(buffers[i], 1, buffer_sizes[i], outfile);
		if (buffer_sizes[i] != nchars) {
			std::cerr << "Expect to write " << buffer_sizes[i] << " but actually wrote " << nchars << std::endl;
			exit(1);
		}
		buffer_sizes[i] = 0;
	}
}

int cmddist_print(FILE *outfile, const Entity &query, const Entity &target, double mutdist, size_t intersize,
		size_t sketchsize64, double mthres, double pthres, size_t raw_intersize, size_t raw_unionsize) {
	if (mutdist > mthres) { return 1; }
	double andprob = query.matchprob * target.matchprob;
	double p = andprob / (query.matchprob + target.matchprob - andprob);
	double pvalue = bhmath_pbinom_tail(intersize, sketchsize64 * NBITS(uint64_t), p);
	if (pvalue > pthres) { return 2; }
	return fprintf(outfile, "%s\t%s\t%f\t%f\t%lu/%lu\n", query.name.c_str(), target.name.c_str(), mutdist, pvalue, raw_intersize, raw_unionsize);
}

const char *ordinal_num_to_suffix(const size_t n) {
	if (n % 10 > 3 || n % 10 == 0 || (10 <= n % 100 && n % 100 <= 20)) {
		return "th";
	}
	if (n % 10 == 1) { return "st"; }
	if (n % 10 == 2) { return "nd"; }
	if (n % 10 == 3) { return "rd"; }
	abort();
}

inline void printinfo(const size_t i2max, const size_t j2max, const size_t i2gridmax, const size_t j2gridmax,
		const clock_t &begclock, const time_t &begtime) {
	if ((ispowof2(i2max) || i2max == i2gridmax) && (ispowof2(j2max) || j2max == j2gridmax)) {
		time_t endtime;
		time(&endtime);
		std::cerr << "Processed up to the cached chunk at (" << i2max << "," << j2max << ") in "
		          << (clock() - begclock) / CLOCKS_PER_SEC << " CPU seconds and "
		          << difftime(endtime, begtime) << " wall-clock seconds." << std::endl;
	}
}

#if DIST_IS_TEMPLATED
template<bool tCLUSTER, bool tNNEIGHBORS>
void cmddist(
#else
void cmddist(bool tCLUSTER, bool tNNEIGHBORS,
#endif
		const std::vector<Entity> &entities1, const std::vector<Entity> &entities2,
		const Sketchargs &args1, const Distargs &args2) {

	size_t nthreads = args2.nthreads;
	if (0 == nthreads) {
#if defined(_OPENMP)
		nthreads = omp_get_max_threads();
#else
		nthreads = 1;
#endif
	}

	FILE *outfile;
	if ("-" == args2.outfname) {
		outfile = stdout;
	} else {
		outfile = fopen(args2.outfname.c_str(), "w");
		if (NULL == outfile) {
			std::cerr << "Unable to open the file " << args2.outfname << " for writing.\n";
			exit(-1);
		}
	}

	std::vector<double> intersize_to_mutdist;
	if (1 == args2.model){
		intersize_to_mutdist_init_poisson(intersize_to_mutdist, args1.sketchsize64, args1.kmerlen);
	} else if (2 == args2.model) {
		intersize_to_mutdist_init_binomial(intersize_to_mutdist, args1.sketchsize64, args1.kmerlen);

	}
	std::cerr << "Max number of genome comparions per cached data: " << ISIZE << "x" << CACHE_SIZE << std::endl;

	size_t buffer_sizes[ISIZE];
	size_t buffer_maxsizes[ISIZE];
	char *buffers[ISIZE];

	for (size_t i = 0; i < ISIZE; i++) {
		buffers[i] = (char*)xmalloc(LINELEN_MAX);
		buffer_sizes[i] = 0;
		buffer_maxsizes[i] = LINELEN_MAX;
	}

	clock_t begclock = clock();
	time_t begtime;
	time(&begtime);

	if (!tNNEIGHBORS) {

		for (size_t i2 = 0; i2 < entities1.size(); i2 += ISIZE) {
			size_t i2max = MIN(i2 + ISIZE, entities1.size());
			size_t j2start = (tCLUSTER ? (i2+0) : 0);
			for (size_t j2 = j2start; j2 < entities2.size(); j2 += CACHE_SIZE) {
					size_t j2max = MIN(j2 + CACHE_SIZE, entities2.size());

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
#endif
				for (size_t i = i2; i < i2max; i++) {
					for (size_t j = (tCLUSTER ? MAX(i+1, j2) : j2); j < j2max; j++) {
						double mutdist, pvalue;
						size_t raw_intersize, raw_unionsize;
						int filter = cmddist_filter(mutdist, pvalue,
								raw_intersize, raw_unionsize,
								entities1[i], entities2[j],
								args1, args2,
								intersize_to_mutdist);
						if (filter) { continue; }
						size_t iter = i - i2;
						assert(iter < ISIZE);
						cmddist_buffer(buffers[iter], buffer_sizes[iter], buffer_maxsizes[iter],
								entities1[i], entities2[j], mutdist, pvalue, raw_intersize, raw_unionsize);
					}
				}
				printinfo(i2max, j2max, entities2.size(), entities2.size(), begclock, begtime);
				cmddist_write(outfile, buffers, buffer_sizes, i2max - i2);
			}
		}

	} else {
		typedef std::tuple<double, double, size_t, size_t, size_t, size_t> hit_t;
		std::array<std::priority_queue<hit_t>, ISIZE> tophits_arr;
		for (size_t i2 = 0; i2 < entities1.size(); i2 += ISIZE) {
			size_t i2max = MIN(i2 + ISIZE, entities1.size());
			for (size_t j2 = 0; j2 < entities2.size(); j2 += CACHE_SIZE) {
				size_t j2max = MIN(j2 + CACHE_SIZE, entities2.size());

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
#endif
				for (size_t i = i2; i < i2max; i++) {
					for (size_t j = j2; j < j2max; j++) {
						if (i == j) { continue; }
						double mutdist, pvalue;
						size_t raw_intersize, raw_unionsize;
						int filter = cmddist_filter(mutdist, pvalue,
								raw_intersize, raw_unionsize,
								entities1[i], entities2[j],
								args1, args2,
								intersize_to_mutdist);
						if (filter) { continue; }
						size_t iter = i - i2;
						assert(iter < ISIZE);
						hit_t hit = std::make_tuple(mutdist, pvalue, raw_intersize, raw_unionsize, i, j);
						auto &tophits = tophits_arr[iter];
						if (tophits.size() == 0 || hit < tophits.top()) {
							tophits.push(hit);
							if (tophits.size() == args2.nneighbors) {
								tophits.pop();
							}
						}
					}
				}
				printinfo(i2max, j2max, entities1.size(), entities2.size(), begclock, begtime);
			}

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
#endif
			for (size_t i = i2; i < i2max; i++) {
				size_t iter = i - i2;
				auto &tophits = tophits_arr[iter];
				while (tophits.size() != 0) {
					auto hit = tophits.top();
					tophits.pop();
					auto mutdist = std::get<0>(hit);
					auto pvalue = std::get<1>(hit);
					auto raw_intersize = std::get<2>(hit);
					auto raw_unionsize = std::get<3>(hit);
					auto qid = std::get<4>(hit);
					auto tid = std::get<5>(hit);
					cmddist_buffer(buffers[iter], buffer_sizes[iter], buffer_maxsizes[iter],
							entities1[qid], entities2[tid], mutdist, pvalue, raw_intersize, raw_unionsize);
				}
			}
			cmddist_write(outfile, buffers, buffer_sizes, i2max - i2);
		}
	}

	for (size_t i = 0; i < ISIZE; i++) {
		free(buffers[i]);
	}
	if (NULL != outfile) { fclose(outfile); }
}

int main(int argc, char **argv) {
	time_t begtime, endtime;
	time(&begtime);

	std::string commitsuffix = "";
	if (std::string("()") == (STR((GIT_DIFF_SHORTSTAT)))) {
		commitsuffix = "-clean";
	} else {
		commitsuffix = std::string("-dirty ") + (STR((GIT_DIFF_SHORTSTAT)));
	}

	if (argc < 2 || !strcmp("--help", argv[1])) { allusage(argc, argv); }
	if (!strcmp("--version", argv[1])) {
		std::cerr << "version " << VERSION << " commit " << (STR(GIT_COMMIT_HASH)) << commitsuffix << std::endl;
		exit(-1);
	}
	std::cerr << "Running " << argv[0] << " commit " << (STR(GIT_COMMIT_HASH)) << commitsuffix << std::endl;


	if (!strcmp("exact", argv[1])) {
		Sketchargs args;
		args.parse(argc, argv);
		std::vector<Kmerset> kmersetVec(args.infnames.size(), Kmerset());
		for (size_t i = 0; i < kmersetVec.size(); i++) {
			kmerset_init(kmersetVec[i], args.infnames[i].c_str(), args.kmerlen, args.iscasepreserved, args.isstrandpreserved);
		}
		for (size_t i = 0; i < kmersetVec.size(); i++) {
			for (size_t j = i + 1; j < kmersetVec.size(); j++) {
				kmerset_dist(kmersetVec[i], kmersetVec[j], args.kmerlen);
			}
		}
	} else if (!strcmp("sketch", argv[1])) {
		Sketchargs args;
		args.parse(argc, argv);

		std::vector<std::vector<std::pair<size_t, size_t>>> fid_to_entityid_count_list;
		std::vector<std::string> fid_to_fname;
		std::vector<std::string> entityid_to_entityname;
		assert (0 < args.infnames.size());
		// parse filecontents
		parse_metaf(fid_to_entityid_count_list, fid_to_fname, entityid_to_entityname, args.infnames);
		assert (fid_to_entityid_count_list.size() == fid_to_fname.size());
		std::vector<Entity> entities(entityid_to_entityname.size(), Entity());


		clock_t para_begclock = clock();
		time_t para_begtime, para_endtime;
		time(&para_begtime);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1) num_threads((args.nthreads ? args.nthreads : omp_get_max_threads()))
#endif
		for (size_t i = 0; i < fid_to_fname.size(); i++) {
			entities_init(entities, fid_to_fname[i], args, fid_to_entityid_count_list[i], entityid_to_entityname);
			// entity_init(entities[i], args.infnames[i], args);
			if (0 == (i & (i + 1))) {
				time(&para_endtime);
				std::cerr << "Sketching the first " <<  i + 1 << " files consumed "
				          << (clock() - para_begclock) / CLOCKS_PER_SEC << " CPU seconds and "
				          << difftime(para_endtime, para_begtime) << " real seconds" << std::endl;
				          // << "The last sequence is converted into " << entities[i].usigs.size() << " 64-bit integers.\n";
			}
		}
		save_entities(entities, args.outfname);
		time(&endtime);
		std::string systime_began2 = ctime(&begtime);
		std::string systime_ended2 = ctime(&endtime);
		args.write(systime_began2, systime_ended2, clock() / CLOCKS_PER_SEC);
	} else if (!strcmp("dist", argv[1])) {
		Distargs args;
		args.parse(argc, argv);

		Sketchargs args1;
		std::vector<Entity> entities1;
		args1.read(args.infnames[0]); // args1.usage(argc, argv);
		load_entities(entities1, args.infnames[0]);

		std::vector<Entity> entities2;
		if (args.infnames.size() == 2) {
			Sketchargs args2;
			args2.read(args.infnames[1]);
			std::string arg;
			verifycompatible(args1, args2, args.infnames[0], args.infnames[1]);
			load_entities(entities2, args.infnames[1]);
			//for (auto e : entities1) { print_entity(e); }
			//for (auto e : entities2) { print_entity(e); }
		}

		const bool tCLUSTER = (1 == args.infnames.size());
		const bool tNNEIGHBORS = (0 < args.nneighbors);
#if DIST_IS_TEMPLATED
		if (tCLUSTER && tNNEIGHBORS) {
			cmddist<true, true>(entities1, entities1, args1, args);
		} else if (!tCLUSTER && tNNEIGHBORS) {
			cmddist<false, true>(entities1, entities2, args1, args);
		} else if (tCLUSTER && !tNNEIGHBORS) {
			cmddist<true, false>(entities1, entities1, args1, args);
		} else {
			cmddist<false, false>(entities1, entities2, args1, args);
		}
#else
		cmddist(tCLUSTER, tNNEIGHBORS, entities1, (tCLUSTER ? entities1 : entities2), args1, args);
#endif
	} else {
		std::cerr << "Unrecognized command: " << argv[1] <<  "\n";
		allusage(argc, argv);
	}

	time(&endtime);
	std::cerr << "Program ran for " << clock() / CLOCKS_PER_SEC << " CPU seconds and " << difftime(endtime, begtime) << " real seconds." << std::endl;
}

