[![CodeQL](https://github.com/jermp/lphash/actions/workflows/codeql.yml/badge.svg)](https://github.com/jermp/lphash/actions/workflows/codeql.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7239205.svg)](https://doi.org/10.5281/zenodo.7239205)

LPHash
======

LPHash is a *minimal perfect hash function* (or MPHF) designed for k-mer sets where overlaps of k-1 bases between consecutive k-mers are exploited to:

1. Reduce the space of the data structure. For example, with the right combination of parameters, LPHash is able to achive < 0.9 bits/k-mer **in practice** whereas the known *theoretical* lower-bound of a classic MPHF is 1.44 bits/k-mer and practical constructions take 2-3 bits/k-mer.

2. Boost the evaluation speed of queries performed for consecutive k-mers.

3. Preserve the *locality* of k-mers: consecutive k-mers are likely to receive consecutive hash codes.

The data structure and its construction/query algorithms are described in the paper *Locality-Preserving Minimal Perfect Hashing of k-mers* [1].
A pre-print is available here: [https://arxiv.org/abs/2210.13097](https://arxiv.org/abs/2210.13097).

Datasets are available on Zenodo here: [https://zenodo.org/record/7239205](https://zenodo.org/record/7239205).

#### Table of contents
* [Compiling the Code](#compiling-the-code)
* [Dependencies](#dependencies)
* [Tools](#tools)
* [Build a Function](#build-a-function)
* [Example](#Example)
* [Input Files](#input-files)
* [Authors](#authors)
* [References](#references)

Compiling the Code
------------------

The code is tested on Linux with `gcc` and on Mac OS with `clang` (Intel and ARM processors supported, like Apple M1).
To build the code, [`CMake`](https://cmake.org/) is required.

Clone the repository with

    git clone --recursive https://github.com/jermp/lphash.git

If you have cloned the repository **without** `--recursive`, be sure you pull the dependencies with the following command before compiling:

    git submodule update --init --recursive

K-mers can have maximal lengths of both 32 or 64, depending on their definition found in the file `include/compile_constants.tpd` (either `uint64_t` or `__uint128_t`). Although larger 128-bit k-mers also work for k <= 32, we recommend to use the right type whenever possible.

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following:

    mkdir build
    cd build
    cmake ..
    make -j

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D LPHASH_USE_SANITIZERS=On
    make -j

Dependencies
------------

The repository has minimal dependencies: it only uses the [PTHash](https://github.com/jermp/pthash) library (for minimal perfect hashing [2,3]), and `zlib` to read gzip-compressed streams.

To automatically pull the PTHash dependency, just clone the repo with
`--recursive` as explained in [Compiling the Code](#compiling-the-code).

If you do not have `zlib` installed, you can do

    sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

    brew install zlib

if you have a Mac (and Homebrew installed: https://brew.sh/).


Tools
-----

There is one main executable called `lphash` after the compilation, which can be used to run a tool.
Run `./lphash` as follows to see a list of available tools.

	LP-Hash: (L)ocality (P)reserving Minimal Perfect (Hash)ing of k-mers

	Usage: ./lphash <tool> ...

	Available tools:
	  build-p      build a partitioned LP-MPHF
	  build-u      build an unpartitioned LP-MPHF
	  query-p      query a partitioned LP-MPHF
	  query-u      query an unpartitioned LP-MPHF


Build a Function
----------------

The driver program called `lphash` can be used to build and query Locality-Preserving MPHFs.
The subcommand `build-p` uses the partitioning strategy, while `build-u` does not (check out the paper for more details).
In the following, we will focus on `build-p` alone since it is the most memory efficient between the two.

From within the directory where the code was compiled (see the section [Compiling the Code](#compiling-the-code)), run the command:

    ./lphash build-p --help

to show the usage of the subcommand (reported below for convenience).

	Usage: build-p [-h,--help] [-i input_filename] [-k k] [-m m] [-s seed] [-t threads] [-o output_filename] [-d tmp_dirname] [-c c] [--max-memory max-memory] [--check] [--verbose]
	
	 [-i input_filename]
		REQUIRED: Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:
		- without duplicate nor invalid kmers
		- one DNA sequence per line.
		For example, it could be the de Bruijn graph topology output by BCALM.
	
	 [-k k]
		REQUIRED: K-mer length (must be <= 63).
	
	 [-m m]
		REQUIRED: Minimizer length (must be < k and <= 32).
	
	 [-s seed]
		Seed for minimizer computation (default is 42).
	
	 [-t threads]
		Number of threads for pthash (default is 1).
	
	 [-o output_filename]
		Output file name where the data structure will be serialized (no files generated by default).
	
	 [-d tmp_dirname]
		Temporary directory used for construction in external memory (default is directory '.').
	
	 [-c c]
		A (floating point) constant that trades construction speed for space effectiveness of minimal perfect hashing. 
		A reasonable value lies between 3.0 and 10.0 (default is 3.00).
	
	 [--max-memory max-memory]
		Maximum internal memory [GB] for building (8GB by default). Use external memory if needed.
	
	 [--check]
		Check correctness after construction (disabled by default).
	
	 [--verbose]
		Verbose output during construction (disabled by default).
	
	 [-h,--help]
		Print this help text and silently exits.

### Output format

`./lphash build-p` prints useful statistics about the constructed MPHF on stdout as single comma-separated values whose fields are:

1. The input filename used to build the MPHF,
2. k-mer length k,
3. minimizer length m,
4. the fraction of non-unique minimizer,
5. the theoretical epsilon (number of super-k-mers / number of k-mers),
6. the true epsilon,
7. alpha (fragmentation factor = number of contigs / number of k-mers),
8. the space of the MPHF in bits/k-mer.

Example
-------

For the examples, we are going to use some collections of *stitched unitigs* from the directory `data/unitigs_stitched`.
These collections were built for k = 31, 47 and 63, so MPHFs should be built with the right k to ensure correctness.

In the section [Input Files](#input-files), we explain how such collections of stitched unitigs can be obtained from raw FASTA files.

This example

    ./lphash build-p -i ../data/unitigs_stitched/se.ust.k31.fa.gz -k 31 -m 15 --check -o se_k31_m15.lph --verbose

builds a dictionary for the k-mers read from the file `data/unitigs_stitched/se.ust.k31.fa.gz`, with k = 31 and m = 15.
It also check the correctness of the MPHF (`--check` option, for both minimality and collisions), and serializes it on disk to the file `se_k31_m15.lph`.


Another example is

	./lphash build-p -i ../data/unitigs_stitched/se.ust.k63.fa.gz -k 63 -m 17 --check -o se_k63_m17.lph --verbose

which uses k = 63 and m = 17.

Similarly to `build-p`, the `query-p` command (or `query-u` if the MPHF was built by `build-u`) outputs single CSV lines on stdout.
Queries can be performed by recomputing each hash value from scratch (random query) of by re-using the fact that successive k-mers overlap by k-1 bases (streaming).
The fields are:
input file,lphash file,total k-mers,barebone streaming time,barebone random time,full streaming time,full random time
1. Input (fasta/fastq) file,
2. the queried MPHF file,
3. the total number of k-mers in the input,
4. average streaming query time in ns/k-mer (without saving the hash values),
5. average random query time in ns/k-mer (without saving the hash values),
6. average streaming query time in ns/k-mer (by saving the hash values inside an in-memory vector),
7. average random query time in ns/k-mer (by saving the hash values inside an in-memory vector).

### Example

    ./lphash query-p -i se_k31_m15.lph -q ../data/queries/salmonella_enterica.fasta.gz

Input Files
-----------

LPHash is meant to index k-mers from collections that do not contain duplicates.
These collections can be obtained, for example, by extracting the maximal unitigs of a de Bruijn graph.

To do so, we can use the tool [BCALM2](https://github.com/GATB/bcalm) [4].
This tool builds a compacted de Bruijn graph and outputs its maximal unitigs.
From the output of BCALM2, we can then *stitch* (i.e., glue) some unitigs to reduce the number of nucleotides.
The stitiching process is carried out using the [UST](https://github.com/jermp/UST) tool.

Below we provide a complete example (assuming both BCALM2 and UST are installed correctly) that downloads the Human (GRCh38) Chromosome 13 and extracts the maximal stitiched unitigs for k = 31.

    mkdir DNA_datasets
    wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -O DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
    ~/bcalm/build/bcalm -in ~/DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 8
    ~/UST/ust -k 31 -i ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa
    gzip Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa.ust.fa
    rm ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa


### Zenodo
We uploaded a copy of all the datasets used in the experiments of the paper, **processed** as we described above, on Zenodo: [https://zenodo.org/record/7239205](https://zenodo.org/record/7239205).

Authors
-------

- [Yoshihiro Shibuya](https://github.com/yhhshb) - <yoshihiro.shibuya@esiee.fr>

- [Giulio Ermanno Pibiri](https://jermp.github.io) - <giulioermanno.pibiri@unive.it>

- [Antoine Limasset](https://github.com/Malfoy) - <antoine.limasset@univ-lille.fr>

References
-----
* [1] Giulio Ermanno Pibiri, Yoshihiro Shibuya, and Antoine Limasset. Locality-Preserving Minimal Perfect Hashing of k-mers. ArXiv. [https://arxiv.org/abs/2210.13097](https://arxiv.org/abs/2210.13097).
* [2] Giulio Ermanno Pibiri and Roberto Trani. [PTHash: Revisiting FCH Minimal Perfect Hashing](https://dl.acm.org/doi/10.1145/3404835.3462849). In The 44th International ACM SIGIR Conference on Research and Development in Information Retrieval, pages 1339–1348, 2021.
* [3] Giulio Ermanno Pibiri and Roberto Trani. [Parallel and external-memory construction of minimal perfect hash functions with PTHash](https://arxiv.org/abs/2106.02350). CoRR, abs/2106.02350, 2021.
* [4] Rayan Chikhi, Antoine Limasset, and Paul Medvedev. [Compacting de Bruijn graphs from sequencing data quickly and in low memory](https://academic.oup.com/bioinformatics/article/32/12/i201/2289008?login=false). Bioinformatics, 32(12):i201–i208, 2016.