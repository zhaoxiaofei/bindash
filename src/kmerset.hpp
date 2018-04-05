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

#ifndef KMERSET_HPP
#define KMERSET_HPP
#include "bhmath.hpp"
#include "cbuffer.hpp"

#include <set>
#include <string>
#include <string.h>

class Kmerset {
public:
	std::string name;
	std::set<std::string> kmers;
	double matchprob;
};

void revcomplement(std::string &rc, const std::string &fw) {	
	for (size_t i = 0; i < fw.size(); i++) {		
		rc[i] = RCMAP[(size_t)fw[fw.size()- 1-i]];
	}
}

void kmerset_init(Kmerset &kmerset, const char *fname, size_t kmerlen, bool iscasepreserved, bool isstrandpreserved) {
	kmerset.name = fname;
	XFILE file = XOPEN(fname, "r");
	assert(file != NULL);
	CBuf cbuf(kmerlen, iscasepreserved);
	while (1) {
		cbuf.eatnext(file);
		if (XEOF(file)) { break; }
		if (kmerlen <= cbuf.slen) {
			std::string kmer = cbuf.tostring();
			if (!isstrandpreserved) {			
				std::string kmer2(kmer.size(), '\0');
				revcomplement(kmer2, kmer);
				if (kmer2 < kmer) {
					kmer = kmer2;
				}
			}
			kmerset.kmers.insert(kmer);
		}
	}
	XCLOSE(file);
	double entropy = bhmath_calc_entropy(cbuf.chfreqs, 256);
	kmerset.matchprob = bhmath_matchprob(kmerlen, entropy, kmerset.kmers.size());
}

template <typename T>
size_t set_calc_intersize(const std::set<T> &setA, const std::set<T> &setB) {
	if (setA.size() > setB.size()) { return set_calc_intersize(setB, setA); }
	size_t intersize = 0;
	for (T elemA : setA) {
		if (setB.find(elemA) != setB.end()) {
			intersize++;
		}
	}
	return intersize;
}

void kmerset_dist(const Kmerset &ks1, const Kmerset &ks2, size_t kmerlen) {
	size_t intersize = set_calc_intersize(ks1.kmers, ks2.kmers);
	size_t unionsize = (ks1.kmers.size() + ks2.kmers.size() - intersize);
	double jaccard_idx = (double)intersize / (double)unionsize;
	double dice = bhmath_jaccard_to_dice(jaccard_idx);
	double mutdist = -1.0 / kmerlen * log(dice);
	double andprob = ks1.matchprob * ks2.matchprob;
	double p = andprob / (ks1.matchprob + ks2.matchprob - andprob);
	//double pvalue = bhmath_pbinom(sketchsize64 * NBITS(uint64_t), p, intersize);
	size_t n = MIN(ks1.kmers.size(), ks2.kmers.size());
	double sd = sqrt(n * p * (1-p));
	printf("%s\t%s\t%.4f\terf(%.4f)\t%lu\t%lu\n",
		   ks1.name.c_str(), ks2.name.c_str(),
		   mutdist, (jaccard_idx - p) * (double)n / sd, intersize, unionsize);
}

#endif
