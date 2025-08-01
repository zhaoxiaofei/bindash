[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bindash/README.html)
![](https://anaconda.org/bioconda/bindash/badges/license.svg)
![](https://anaconda.org/bioconda/bindash/badges/version.svg)
![](https://anaconda.org/bioconda/bindash/badges/latest_release_relative_date.svg)
![](https://anaconda.org/bioconda/bindash/badges/platforms.svg)
[![install with conda](https://anaconda.org/bioconda/bindash/badges/downloads.svg)](https://anaconda.org/bioconda/bindash)

<div align="center">
  <img width="40%" src ="BinDash_logo.svg">
</div>

BinDash is a command-line software for comparing genomes (including metagenomes and pangenomes) on a typical personal laptop. BinDash is based on **Bin**wise **D**ensified minh**ash** for estimation of mutation rate between genomes. We implemented ***b-bit one-permutation rolling MinHash with optimal/faster/re-randomized densification***.  It is extremely fast and memory efficient. It can handle sequences consisting of terabytes of input data (gzipped or not, in fasta or fastq format). A Rust implementation can be found [here](https://github.com/jianshu93/bindash-rs).

The basic idea is: the Jaccard index as an accurate proxy of Average Nucleotide Identity(ANI) or mutation rate (1-ANI) according to equation:

$$ANI=1+\frac{1}{k}log\frac{2*J}{1+J}$$

if assuming a Poisson model, and 

$$ANI=(\frac{2*J}{1+J})^{\frac{1}{k}}$$

if assuming a Binomial model. You can specify which model to use via --model option, see below.


If you find BinDash useful, please cite the following papers:

```bash
@article{zhao2019bindash,
  title={BinDash, software for fast genome distance estimation on a typical personal laptop},
  author={Zhao, XiaoFei},
  journal={Bioinformatics},
  volume={35},
  number={4},
  pages={671--673},
  year={2019},
  publisher={Oxford University Press}
}

@article{zhao2024bindash,
  title={Bindash 2.0: new MinHash scheme allows ultra-fast and accurate genome search and comparisons},
  author={Zhao, Jianshu and Zhao, Xiaofei and Pierre-Both, Jean and Konstantinidis, Konstantinos T},
  journal={bioRxiv},
  pages={2024--03},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}

```


# How to install (simple):
Precompiled binaries on modern Linux can be found [here](https://github.com/jianshu93/bindash/releases/tag/v2.1). On MacOS, GNU GCC has to be installed first, we recommend the homebrew install (see below).

```bash
wget https://github.com/jianshu93/bindash/releases/download/v2.1/BinDash_Linux_x86-64_v2.0.tar.gz
tar -xzvf BinDash_Linux_x86-64_v2.0.tar.gz
./bindash --help
```

## If you have conda installed on linux

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bindash/README.html)

```bash
conda install bindash -c bioconda

```
## if you have homebrew insalled on MacOS
```bash
brew tap jianshu93/bindash
brew update
brew install --cc=gcc-13 bindash
```

# How to install from source:
```sh
cd ${PROJECT_ROOT_DIRECTORY}  
mkdir release && cd release
cmake -DCMAKE_BUILD_TYPE=Release ..  
make # For Windows with MSYS Makefiles, the command might be "cd ../ && make" because out-of-source build may or may not be supported on this platform. 
./bindash --help # to see a general help message   
```
For MacOS, the native clang compiler cannot compile without compiling error. It is recommended to install gcc first as follows.

```
# install homebrew if not already done
brew update
brew install gcc@13

cd ${PROJECT_ROOT_DIRECTORY}  
mkdir release && cd release

# Then run cmake as above (a GUI for cmake may also be available for MacOS)
CC="$(brew --prefix)/bin/gcc-13" CXX="$(brew --prefix)/bin/g++-13" cmake -DCMAKE_INSTALL_PREFIX=. ..
make
```
# Dependencies:

 - any C++ compiler supporting the C++11 standard
 - CMake version 2.6 or plus
 - zlib 

# How to run:

The folowing three commands show how to run BinDash:
```sh
bindash sketch --outfname=genomeA.sketch genomeA.fasta
bindash sketch --outfname=genomeB.sketch genomeB.fasta
bindash dist genomeA.sketch genomeB.sketch # print to stdout

### or if you have a list of genomes, one genome per line in a list file. All versus all
ls *.fasta > name.txt
bindash sketch --listfname=name.txt --outfname=genome_sketch
bindash dist --outfname=dist.txt genome_sketch

### query against database genomes
ls *_query.fasta > query.txt
ls *_db.fasta > db.txt
bindash sketch --listfname=query.txt --outfname=genome_query_sketch
bindash sketch --listfname=db.txt --outfname=genome_db_sketch
bindash dist --outfname=dist.txt genome_query_sketch genome_db_sketch
```

The output of "bindash dist" consists of several lines. 
Each line has these five tab-separated fields: 
 - query genome (Q)
 - target genome (T) 
 - mutation distance between Q and T
 - p-value for the mutation distance
 - Jaccard Index between Q and T

According to our experiments and theoretical analysis, sketch size (the skethchsize64 option) shoule be larger than 188 (that is actual sketch size is larger than ~12,000) to be accurate for genome pairs with ANI above 99.5%

# How BinDash works:

By default, BinDash uses the optimally densified MinHash proposed by Shrivastava. This MinHash technique allows for efficient compression and fast comparison. 

Basically, compression of a genome is done as follows.
 1. All k-mers of each genome are selected and then transformed into hash values.
 2. The range of all possible hash values are partitioned into some bins.
 3. The smallest hash value in each bin is selected.
 4. If a bin is empty (i.e., has no hash values), then the smallest hash value from the next non-empty bin is picked (--dens=1). The definition of the next bin is proposed by Shrivastava 2017 (or called optimal densification). Another densification strategy is to map reversely from non-empty bins to empty bins and choose the smallest hash value from non-empty bin to fill the empty bin, as proposed by Mai et.al. 2020 (--dens=2). This new densification strategy has a worse case O(k*log(k)) complexity  while Shrivastava 2017 has a O(k^2) worse case complexity. The last densification strategy is called Re-randomized densification (--dens=3), which rans an additional MinHash step within each previously empty bin after being filled via the optimal densification, see Li et.al. 2019. 
 5. The lowest (i.e., least significant) b-bits of each hash value are picked to form the signature of the genome.
 6. Two genomes are compared by simply performing b XOR opeations for b bit positions, followed by (b-1) AND operations for these b bit positions. 
 7. By default, BinDash automatically picks the fastest POPCOUNT algorithm for the given array or sketch size based on your machine instructions (e.g. AVX2 or AVX512).


# Limitations:

- The k-mer size has to increase as the corresponding genome size increases. A natural solution is to use k-mer frequency instead of k-mer presence/abasence. However, this natural solution relies on a probabilistic model for k-mer frequency. Such model may be heavily genome-dependent.
- Some large-scale genome rearrangement (e.g., chromosome duplication) brings very little change to Jaccard Index and mutation distance. Weighted Jaccard index estimation via BagMinHash or DartMinHash can be very useful chromosome duplication (no positional information) while Order MinHash (assume linear, non-fragmented genomes) can be used for kmer positional information.

# Other implementations

Rust implementation can be found [here](https://github.com/jianshu93/bindash-rs)

# Additional Information:

All suggestions, comments, and feature requests are welcome.

Author: XiaoFei Zhao (cndfeifei AT hotmail DOT com)  and Jianshu Zhao (jianshuzhao AT yahoo DOT com)
License: Apache 2.0

# Reference

Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal (2014) 57 (11): 1624-1638. http://arxiv.org/abs/1202.4961

Daniel Lemire, Owen Kaser: Recursive n-gram hashing is pairwise independent, at best, Computer Speech & Language, Volume 24, Issue 4, October 2010, Pages 698-710 http://arxiv.org/abs/0705.4676

Daniel Lemire, The universality of iterated hashing over variable-length strings, Discrete Applied Mathematics 160 (4-5), 2012. http://arxiv.org/abs/1008.1715

Muła, Wojciech, Nathan Kurz, and Daniel Lemire. "Faster population counts using AVX2 instructions." The Computer Journal 61, no. 1 (2018): 111-120. https://academic.oup.com/comjnl/article-abstract/61/1/111/3852071 

Anshumali Shrivastava, Optimal Densification for Fast and Accurate Minwise Hashing. Proceedings of the 34th International Conference on Machine Learning, PMLR 70:3154-3163, 2017. http://proceedings.mlr.press/v70/shrivastava17a.html 

Tung Mai et.al., On Densification for Minwise Hashing, Uncertainty in Artificial Intelligence, 2020. http://proceedings.mlr.press/v115/mai20a.html 

Ping Li et.al., Re-randomized Densification for One Permutation Hashing and Bin-wise Consistent Weighted Sampling, 33rd Conference on Neural Information Processing Systems (NIPS), 2019. https://proceedings.neurips.cc/paper/2019/hash/9f067d8d6df2d4b8c64fb4c084d6c208-Abstract.html

XiaoFei Zhao; BinDash, software for fast genome distance estimation on a typical personal laptop. Bioinformatics, 2018. bty651, https://doi.org/10.1093/bioinformatics/bty651

Zhao, J., Zhao, X., Pierre-Both, J. and Konstantinidis, K.T., 2024. Bindash 2.0: new MinHash scheme allows ultra-fast and accurate genome search and comparisons. bioRxiv, pp.2024-03. https://www.biorxiv.org/content/10.1101/2024.03.13.584875v1.abstract
