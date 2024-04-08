#include <vector>
#include <unordered_set>
#include <random>
#include <cstdint>
#include <iostream>
#include "../murmur3.h" // Include MurmurHash3 library
#include "../xxhash_header_only.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Utility function to hash an integer using MurmurHash3, for mapping non-empty bins to empty ones
uint32_t murmur_hash(uint32_t value, uint32_t seed) {
    uint32_t hash_some[1];
    MurmurHash3_x86_32(&value, sizeof(value), seed, hash_some);
    return hash_some[0];
}


/* This densitfication strategy is found at
@article{Mai2020Densification,
  title={On Densification for Minwise Hashing},
  author={Mai, Tung},
  journal={Uncertainty in Artificial Intelligence, http://proceedings.mlr.press/v115/mai20a.html},
  year={2020}
}
*/
int revopt_densify(std::vector<uint64_t>& signs) {

	uint64_t minval = UINT64_MAX;
	uint64_t maxval = 0;
	for (auto sign : signs) { 
		minval = MIN(minval, sign);
		maxval = MAX(maxval, sign);
	}
	if (UINT64_MAX != maxval) { return 0; }
	if (UINT64_MAX == minval) { return -1; }

    size_t size = signs.size();
    std::unordered_set<size_t> E; // Set of empty bins

    // Initial identification of empty bins, we use UINT64_MAX as unassigned bins
    for (size_t i = 0; i < size; ++i) {
        if (signs[i] == UINT64_MAX) {
            E.insert(i);
        }
    }
    // Main densification loop
    size_t alpha = 0;
    while (!E.empty()) {
        for (size_t j = 0; j < size; ++j) {
            if (signs[j] != UINT64_MAX) {
                size_t i = murmur_hash(j, alpha) % size;
                if (E.find(i) != E.end()) {
                    signs[i] = signs[j]; // Directly update the bin in the input vector
                    E.erase(i);
                    if (E.empty()) break;
                }
            }
        }
        alpha++;
    }
    return 1; // Return a status code indicating success
}

uint64_t univhash2(uint64_t s, uint64_t t) {
    uint64_t x = (1009) * s + (1000000003ULL) * t;
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
int opt_densify(std::vector<uint64_t>& signs) {
    uint64_t minval = UINT64_MAX;
    uint64_t maxval = 0;
    for (auto sign : signs) { 
        minval = MIN(minval, sign);
        maxval = MAX(maxval, sign);
    }
    if (UINT64_MAX != maxval) { return 0; }
    if (UINT64_MAX == minval) { return -1; }
    for (uint64_t i = 0; i < signs.size(); i++) {
        if (signs[i] == UINT64_MAX) {
            uint64_t j = i;
            uint64_t nattempts = 0;
            while (UINT64_MAX == signs[j]) {
                j = univhash2(i, nattempts) % signs.size();
                nattempts++;
            }
            signs[i] = signs[j];
        }
    }
    return 1;
}

/* This densitfication strategy is found at
@article{Li2019re-randomized,
  title={Re-randomized Densification for One Permutation Hashing and Bin-wise Consistent Weighted Sampling},
  author={Li, Ping},
  journal={Neural Information Processing Systems, 2019},
  year={2019}
}
*/
int rerand_densify(std::vector<uint64_t>& signs) {
    uint64_t minval = UINT64_MAX;
    uint64_t maxval = 0;
    for (auto sign : signs) {
        minval = MIN(minval, sign);
        maxval = MAX(maxval, sign);
    }
    if (UINT64_MAX != maxval) { return 0; }
    if (UINT64_MAX == minval) { return -1; }

    for (uint64_t i = 0; i < signs.size(); i++) {
        if (signs[i] == UINT64_MAX) {
            uint64_t j = i;
            uint64_t nattempts = 0;
            while (UINT64_MAX == signs[j]) {
                j = univhash2(i, nattempts) % signs.size();
                nattempts++;
            }
            // Using xxHash for the MinHash step before filling the empty bin. One hash function to approximate the permutation
            uint64_t seed = i;  // Using i as the seed for variability
            signs[i] = XXH64(&signs[j], sizeof(signs[j]), seed);
        }
    }
    return 1;
}


// generate some large sparse vectors
std::vector<uint64_t> generateSparseVector(size_t size, double nonEmptyRatio) {
    std::vector<uint64_t> vec(size, UINT64_MAX);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (size_t i = 0; i < size; ++i) {
        if (distribution(generator) < nonEmptyRatio) {
            vec[i] = i; // Assigning a non-empty value (can be any value other than UINT64_MAX)
        }
    }

    return vec;
}

int main() {
    size_t vectorSize = 100;
    double nonEmptyRatio = 0.85; // 85% of the elements will be non-empty

    std::vector<uint64_t> sparse_vector = generateSparseVector(vectorSize, nonEmptyRatio);
    std::vector<uint64_t> sparse_vector_copy = sparse_vector; // Make a copy for the second densification function
    std::vector<uint64_t> sparse_vector_copy2 = sparse_vector; // Make a copy for the second densification function


    int status_opt = opt_densify(sparse_vector);
    int status_revopt = revopt_densify(sparse_vector_copy);
    int status_rerand = rerand_densify(sparse_vector_copy2);
    // Print densified vector from opt_densify
    if (status_opt != -1) {
        std::cout << "Opt Densified Vector: ";
        for (uint64_t value : sparse_vector) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "Opt Densified Vector: Empty or All bins were empty" << std::endl;
    }

    // Print densified vector from revopt_densify
    if (status_revopt == 1) {
        std::cout << "RevOpt Densified Vector: ";
        for (uint64_t value : sparse_vector_copy) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "RevOpt Densified Vector: Operation failed" << std::endl;
    }
    // Print re-randomized densified vector from revopt_densify
    if (status_rerand == 1) {
        std::cout << "Re-randomized Densified Vector: ";
        for (uint64_t value : sparse_vector_copy2) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "Rerand Densified Vector: Operation failed" << std::endl;
    }
    return 0;
}
