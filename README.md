# BinHash
Fast and precise comparison of genomes and metagenomes (in the order of terabytes) on a typical personal laptop

BinHash is a command-line software for comparing genomes (including metagenomes and pangenomes) on a typical personal laptop. 
BinHash is extremely fast and memory efficient.
BinHash can handle sequences consisting of terabytes of input data (gzipped or not, in fasta or fastq format). 
 
Dependencies:
 - any C++ compiler supporting the C++11 standard
 - CMake version 2.6
 - zlib 
 
How to install:

cd ${PROJECT_ROOT_DIRECTORY}  
mkdir release && cd release  
cmake -DCMAKE_BUILD_TYPE=Release ..  
make  
./binhash --help # to see a general help message   


Author: XiaoFei Zhao
License: Apache 2.0
