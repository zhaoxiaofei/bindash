#include "bhmath.hpp"
#include "cbuffer.hpp"
#include "commands.hpp"
#include "genome.hpp"
#include "hashutils.hpp"
#include "kmerset.hpp"

#include "rollinghashcpp/cyclichash.h"

//#include "kseq.h"

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

void entity_init(Entity &entity, const char *fname, const Sketchargs &args) {
	entity.name = fname;
	uint64_t genome_size = 0;

	XFILE file = XOPEN(fname, "r");
	if (NULL == file) { 
		std::cerr << "Unable to open the file " << fname << std::endl; 
		exit(-1);
	}
	
	std::minstd_rand0 g1(args.randseed);
	auto s1 = g1();
	auto s2 = g1();
	CBuf cbuf(args.kmerlen, args.iscasepreserved);
	
	if (-1 == args.minhashtype) { // special perfect hash function for nucleotides
		std::set<uint64_t, std::greater<uint64_t>> signset;
		_DnaPerfectHash hf(args.kmerlen);
		_DnaPerfectHash hfrc(args.kmerlen);
		while (1) {
			cbuf.eatnext(file);
			if (XEOF(file)) { break; }
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
		while (1) {
			cbuf.eatnext(file);
			if (XEOF(file)) { break; }
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
		while (1) {
			cbuf.eatnext(file);
			if (XEOF(file)) { break; }
			hashupdate1(cbuf, signs, hfs, hfrcs, args.isstrandpreserved, args.sketchsize64);
		}
		genome_size = estimate_genome_size1(signs);	
		fillusigs(entity, signs, args.bbits);
	} else if (2 == args.minhashtype) { // multi binvals one hash function
		entity.usigs = std::vector<uint64_t>(args.sketchsize64 * args.bbits, 0);
		std::vector<uint64_t> signs(args.sketchsize64 * NBITS(uint64_t), UINT64_MAX); // carry over
		CyclicHash<uint64_t> hf(args.kmerlen, s1, s2, 64);
		CyclicHash<uint64_t> hfrc(args.kmerlen, s1, s2, 64);
		hashinit2(cbuf, hf, hfrc, args.kmerlen);
		while (1) {
			cbuf.eatnext(file);
			if (XEOF(file)) { break; }
			hashupdate2(cbuf, signs, hf, hfrc, args.isstrandpreserved, args.sketchsize64);
		}
		genome_size = estimate_genome_size2(signs, args.sketchsize64);
		fillusigs(entity, signs, args.bbits);
	}
	double entropy = bhmath_calc_entropy(cbuf.chfreqs, 256);
	entity.matchprob = bhmath_matchprob(args.kmerlen, entropy, genome_size + 1);
	XCLOSE(file);
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

template<bool tCLUSTER, bool tNNEIGHBORS>
void cmddist(const std::vector<Entity> &entities1, const std::vector<Entity> &entities2, 
		const Sketchargs &args1, const Distargs &args) {

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
			outfiles[i] = fopen((args.outfname + "." + std::to_string(i) + ordinal_num_to_suffix(i+1) + args.suffix).c_str(), "w");
			if (NULL == outfiles[i]) {
				std::cerr << "Unable to open the file " << args.outfname << " for writing.\n";
				exit(-1);
			}
		}
	}
	std::vector<double> intersize_to_mutdist;
	intersize_to_mutdist_init(intersize_to_mutdist, args1.sketchsize64, args1.kmerlen);

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
				raw_unionsize = 2 * NBITS(uint64_t) * args1.sketchsize64 - raw_intersize;
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
	}
	for (size_t i = 0; i < nthreads; i++) {
		if (outfiles[i] != stdout) {
			fclose(outfiles[i]);
		}
	}
}


int main(int argc, char **argv) {
	std::cerr << argv[0] << " revision " << (STR(GIT_COMMIT_HASH)) << " " << (STR((GIT_DIFF_SHORTSTAT))) << std::endl;
	auto t = clock();
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
		std::vector<Entity> entities(args.infnames.size(), Entity());

#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1) num_threads((args.nthreads ? args.nthreads : omp_get_max_threads()))
#endif
		for (size_t i = 0; i < entities.size(); i++) {
			entity_init(entities[i], args.infnames[i].c_str(), args);
			if (0 == (i & (i + 1))) {
				std::cerr << "Initialization of the first " <<  i + 1 << " entities consumed " << (clock()-t) / CLOCKS_PER_SEC << " seconds. " 
				          << entities[i].usigs.size() << " kmers were hashes\n";
			}
		}
		save_entities(entities, args.outfname);
		args.write();
		std::cerr << "Hash initialization consumed " << (clock() - t) / CLOCKS_PER_SEC << " seconds" << std::endl;
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
		if (tCLUSTER && tNNEIGHBORS) {
			cmddist<true, true>(entities1, entities1, args1, args);
		} else if (!tCLUSTER && tNNEIGHBORS) {
			cmddist<false, true>(entities1, entities2, args1, args);
		} else if (tCLUSTER && !tNNEIGHBORS) {
			cmddist<true, false>(entities1, entities1, args1, args);
		} else {
			cmddist<false, false>(entities1, entities2, args1, args);
		}
	} else {
		std::cerr << "Unrecognized command: " << argv[1] <<  "\n";
		allusage(argc, argv);
	}
}

