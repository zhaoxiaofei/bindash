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
#include <vector>

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
		int res = densifybin(signs);
		if (res != 0) {
			std::cerr << "Warning: the genome " << entityname << " is densified with flag " << res <<  std::endl;
		}
		genome_size = estimate_genome_size2(signs, args.sketchsize64);
		fillusigs(entity, signs, args.bbits);
	}
	double entropy = bhmath_calc_entropy(cbuf.chfreqs, 256);
	entity.matchprob = bhmath_matchprob(args.kmerlen, entropy, genome_size + 1);
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

int cmddist_print(FILE *outfile, const Entity &query, const Entity &target, double mutdist, size_t intersize,
		size_t sketchsize64, double mthres, double pthres, size_t raw_intersize, size_t raw_unionsize) {
	if (mutdist > mthres) { return 1; }
	double andprob = query.matchprob * target.matchprob;
	double p = andprob / (query.matchprob + target.matchprob - andprob);
	double pvalue = bhmath_pbinom_tail(intersize, sketchsize64 * NBITS(uint64_t), p);
	if (pvalue > pthres) { return 2; }
	return fprintf(outfile, "%s\t%s\t%f\t%f\t%lu/%lu\n", query.name.c_str(), target.name.c_str(), mutdist, pvalue, raw_intersize, raw_unionsize);
#if 0
	char buffer[1024];
	int seqlen = 1000; // snprintf(buffer, 1024, "%s\t%s\t%f\t%f\t%lu/%lu\n", query.name.c_str(), target.name.c_str(), mutdist, pvalue, raw_intersize, raw_unionsize);
	int nchars = MIN(seqlen, 1024);
	size_t ret = -1; // fwrite_unlocked(buffer, 1, nchars, outfile);
	return ret;
#endif
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

/* TODO:
 * Potential improvement to be made: CACHE_AWARE nearest neighbor search (may involve substantial additional coding).
 */
#define CACHE_SIZE (1024*2) // approximated

#define DIST_IS_TEMPLATED 0
#if DIST_IS_TEMPLATED
template<bool tCLUSTER, bool tNNEIGHBORS>
void cmddist(
#else
void cmddist(bool tCLUSTER, bool tNNEIGHBORS,
#endif
		const std::vector<Entity> &entities1, const std::vector<Entity> &entities2, 
		const Sketchargs &args1, const Distargs &args) 
{

	size_t nthreads = args.nthreads;
	if (0 == nthreads) {
#if defined(_OPENMP)
		nthreads = omp_get_max_threads();
#else
		nthreads = 1;
#endif
	}
	std::vector<FILE*> outfiles(nthreads, NULL);
	for (size_t i = 0; i < nthreads; i++) {
		if ("-" == args.outfname) { 
			outfiles[i] = stdout;
		} else {
			outfiles[i] = fopen((args.outfname + "." + std::to_string(i+1) + ordinal_num_to_suffix(i+1) + args.suffix).c_str(), "w");
			if (NULL == outfiles[i]) {
				std::cerr << "Unable to open the file " << args.outfname << " for writing.\n";
				exit(-1);
			}
		}
	}
	std::vector<double> intersize_to_mutdist;
	intersize_to_mutdist_init(intersize_to_mutdist, args1.sketchsize64, args1.kmerlen);
	auto t = clock();

if (!tNNEIGHBORS && CACHE_SIZE > 0 
		// && nthreads > 1
		) {

for (size_t i2 = 0; i2 < entities1.size(); i2 += CACHE_SIZE) {
	size_t i2max = MIN(i2 + CACHE_SIZE, entities1.size());
	for (size_t j2 = 0; j2 < entities2.size(); j2 += CACHE_SIZE) {
		size_t j2max = MIN(j2 + CACHE_SIZE, entities2.size());

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
#endif
		for (size_t i = i2; i < i2max; i++) {
			for (size_t j = (tCLUSTER ? MAX(i+1, j2) : j2); j < j2max; j++) {
				size_t raw_intersize, raw_unionsize, intersize;
				if (args1.minhashtype <= 0) {
					intersize = calc_intersize0(raw_intersize, raw_unionsize, entities1[i], entities2[j], args1.sketchsize64);
				} else {
					intersize = calc_intersize12(entities1[i], entities2[j], args1.sketchsize64, args1.bbits);
					raw_intersize = intersize;
					raw_unionsize = NBITS(uint64_t) * args1.sketchsize64;
				}
				// const size_t interdiff = NBITS(uint64_t) * args1.sketchsize64 - intersize;
				cmddist_print(outfiles[i%nthreads], entities1[i], entities2[j], intersize_to_mutdist[intersize], intersize,
						args1.sketchsize64, args.mthres, args.pthres, raw_intersize, raw_unionsize);
				
			}
		}
		if (0 == (i2 & (i2 + 1)) || 0 == (j2 & (j2+1))) {
			std::cerr << "Processed cache chunk (" << i2 << "," << j2 << ") with size " << CACHE_SIZE << " in " 
			          << (clock() - t) / CLOCKS_PER_SEC << " seconds." << std::endl;
			// if (i2+1 == 1024*8) { exit(0); }
		}

	}
}

} else {

#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1) num_threads(nthreads)
#endif
	for (size_t i = 0; i < entities1.size(); i++) {
		std::priority_queue<std::tuple<size_t, size_t, size_t, Entity>> tophits;
		for (size_t j = ((tCLUSTER && !tNNEIGHBORS) ? (i+1) : 0); j < entities2.size(); j++) {
			size_t raw_intersize, raw_unionsize, intersize;
			if (args1.minhashtype <= 0) {
				intersize = calc_intersize0(raw_intersize, raw_unionsize, entities1[i], entities2[j], args1.sketchsize64);
			} else {
				intersize = calc_intersize12(entities1[i], entities2[j], args1.sketchsize64, args1.bbits);
				raw_intersize = intersize;
				raw_unionsize = NBITS(uint64_t) * args1.sketchsize64;
			}
			const size_t interdiff = NBITS(uint64_t) * args1.sketchsize64 - intersize;
			if (tNNEIGHBORS)	{
				if (tophits.size() == args.nneighbors && std::get<0>(tophits.top()) < interdiff) {
					continue;
				}
				tophits.push(std::make_tuple(interdiff, raw_intersize, raw_unionsize, entities2[j]));
				if (tophits.size() == args.nneighbors + 1) {
					tophits.pop();
				}
			} else {
				cmddist_print(outfiles[i%nthreads], entities1[i], entities2[j], intersize_to_mutdist[intersize], intersize,
						args1.sketchsize64, args.mthres, args.pthres, raw_intersize, raw_unionsize);
			}
		}
		if (tNNEIGHBORS) {
			std::vector<std::tuple<size_t, size_t, size_t, Entity>> targets;
			targets.reserve(tophits.size());			
			while (tophits.size() > 0) {
				targets.push_back(tophits.top());
				tophits.pop();
			}
			for (auto it = targets.rbegin(); it != targets.rend(); it++) {
				size_t intersize = NBITS(uint64_t) * args1.sketchsize64 - std::get<0>(*it);
				cmddist_print(outfiles[i%nthreads], entities1[i], std::get<3>(*it), intersize_to_mutdist[intersize], intersize,
						args1.sketchsize64, args.mthres, args.pthres, std::get<1>(*it), std::get<2>(*it));
			}
		}
		if (0 == (i & (i + 1))) {
			std::cerr << "Processed " << i + 1 << " queries in " << (clock() - t) / CLOCKS_PER_SEC << " seconds." << std::endl;
			// if (i+1 == 1024*8) { exit(0); }
		}
	}
}

	std::cerr << "Distance calculation consumed " << (clock() - t) / CLOCKS_PER_SEC << " seconds." << std::endl;
	
	for (size_t i = 0; i < nthreads; i++) {
		if (outfiles[i] != stdout) {
			fclose(outfiles[i]);
		}
	}
}


int main(int argc, char **argv) {
	std::chrono::system_clock::time_point systime_start = std::chrono::system_clock::now();
	auto t = clock();
	time_t systime_began1, systime_ended1;
	std::cerr << argv[0] << " revision " << (STR(GIT_COMMIT_HASH)) << " " << (STR((GIT_DIFF_SHORTSTAT))) << std::endl;
	if (argc < 2 || !strcmp("--help", argv[1])) { allusage(argc, argv); }
	
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

#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1) num_threads((args.nthreads ? args.nthreads : omp_get_max_threads()))
#endif
		for (size_t i = 0; i < fid_to_fname.size(); i++) {
			entities_init(entities, fid_to_fname[i], args, fid_to_entityid_count_list[i], entityid_to_entityname);
			// entity_init(entities[i], args.infnames[i], args);
			if (0 == (i & (i + 1))) {
				std::cerr << "Sketching the first " <<  i + 1 << " files consumed " << (clock()-t) / CLOCKS_PER_SEC << " seconds. " 
				          << "The last sequence is converted into " << entities[i].usigs.size() << " 64-bit integers.\n";
			}
		}
		save_entities(entities, args.outfname);
		std::chrono::system_clock::time_point systime_end = std::chrono::system_clock::now();
		systime_began1 = std::chrono::system_clock::to_time_t(systime_start);
		systime_ended1 = std::chrono::system_clock::to_time_t(systime_end);
		std::string systime_began2 = ctime(&systime_began1);
		std::string systime_ended2 = ctime(&systime_ended1);	
		args.write(systime_began2, systime_ended2);
		std::cerr << "In total, sketching consumed " << (clock() - t) / CLOCKS_PER_SEC << " seconds" << std::endl;
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
		cmddist(tCLUSTER, tNNEIGHBORS, entities1, entities1, args1, args);
#endif
	} else {
		std::cerr << "Unrecognized command: " << argv[1] <<  "\n";
		allusage(argc, argv);
	}
}

