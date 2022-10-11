LPHash
======

#### Table of contents
* [Compiling the Code](#compiling-the-code)
* [Dependencies](#dependencies)
* [Build a Dictionary](#build-a-dictionary)
* [Examples](#Examples)
* [Input Files](#input-files)
* [Large-Scale Benchmark](#large-scale-benchmark)
* [Author](#author)
* [References](#references)

Compiling the Code
------------------

The code is tested on Linux with `gcc` and on Intel Macs with `clang`.
**Apple Silicon Macs are not yet supported.**
To build the code, [`CMake`](https://cmake.org/) is required.

Clone the repository with

    git clone --recursive https://github.com/jermp/lphash.git

If you have cloned the repository **without** `--recursive`, be sure you pull the dependencies with the following command before compiling:

    git submodule update --init --recursive

K-mers can have maximal lengths of both 32 or 64, depending on their definition found in the file `include/compile_constants.tpd` (either uint64_t or kmer128_t). While larger 128 bit k-mers also work for k <= 32 without problems, we recommend to use the right type whenever possible, for maximum optimization.

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following:

    mkdir build
    cd build
    cmake ..
    make -j

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D SSHASH_USE_SANITIZERS=On
    make -j

Dependencies
------------

The repository has minimal dependencies: it only uses the [PTHash](https://github.com/jermp/pthash) library (for minimal perfect hashing), and `zlib` to read gzip-compressed streams.

To automatically pull the PTHash dependency, just clone the repo with
`--recursive` as explained in [Compiling the Code](#compiling-the-code).

If you do not have `zlib` installed, you can do

    sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

    brew install zlib

if you have a Mac (and Homebrew installed: https://brew.sh/).

Build a Locality-Preserving MPHF
------------------

The driver program called `build` and `build_alt` can be used to build Locality-Preserving MPHFs.
`build` uses the partitioning strategy while `build_alt` does not (check out the paper for more details).
In the following we will focus on `build` alone since it is the most memory efficient between the two.

From within the directory where the code was compiled (see the section [Compiling the Code](#compiling-the-code)), run the command:

    ./build --help

to show the usage of the driver program (reported below for convenience).

    Usage: ./build [-h,--help] input_filename k m [-s seed] [-t threads] [-o output_filename] [-d tmp_dirname] [--check] [--verbose]

    input_filename
        Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:
        - without duplicate nor invalid kmers
        - one DNA sequence per line.
        For example, it could be the de Bruijn graph topology output by BCALM.

    k
        K-mer length (must be <= <K-max>).

    m
        Minimizer length (must be < k).

    [-s seed]
        Seed for minimizer computation (default is 42).

    [-t threads]
        Number of threads for pthash (default is 1).

    [-o output_filename]
        Output file name where the data structure will be serialized.

    [-d tmp_dirname]
        Temporary directory used for construction in external memory. Default is directory '.'.

    [--check]
        Check correctness after construction.

    [--verbose]
        Verbose output during construction.

    [-h,--help]
        Print this help text and silently exits.

### Output format

`build` prints useful statistics about the constructed MPHF on stdout as single comma-separated values whose fields are:

1. The input filename used to build the MPHF,
2. k-mer length k,
3. minimizer length m,
4. the fraction of non-unique minimizer,
5. the theoretical epsilon (number of super-k-mers / number of k-mers),
6. the true epsilon,
7. alpha (fragmentation factor = number of contigs / number of k-mers),
8. the space of the MPHF in bits/k-mer.

Example
--------

For the examples, we are going to use some collections of *stitched unitigs* from the directory `data/unitigs_stitched`.
These collections were built for k = 31, 47 and 63, so MPHFs should be built with the right k to ensure correctness.

In the section [Input Files](#input-files), we explain how such collections of stitched unitigs can be obtained from raw FASTA files.

### Example

    ./build ../data/unitigs_stitched/se.ust.k31.fa.gz 31 15 --check -o -o se_k31_m15.lph

This example builds a dictionary for the k-mers read from the file `data/unitigs_stitched/se.ust.k31.fa.gz`, with k = 31 and m = 15. 
It also check the correctness of the MPHF (`--check` option, for both minimality and collisions), and serializes it on disk to the file `se_k31_m15.lph`.

To run a performance benchmark after construction of the index,
use:

    ./bench salmonella_enterica.index

Query a MPHF
--------

Similarly to `build`, the `query` command (or `query_alt` if the MPHF was built by `build_alt`) outputs single CSV lines on stdout.
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

    ./query se_k31_m15.lph ../data/queries/salmonella_enterica.fasta.gz


<!-- ### Example 3

    ./build ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz 31 13 -l 4 -s 347692 --canonical-parsing -o salmonella_100.canon.index

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (same used in Example 2), with k = 31, m = 13, l = 4, using a seed 347692 for construction (`-s 347692`), and with the canonical parsing modality (option `--canonical-parsing`). The dictionary is serialized on disk to the file `salmonella_100.canon.index`.

The  "canonical" version of the dictionary offers more speed for only a little space increase (for a suitable choice of parameters m and l), especially under low-hit workloads -- when the majority of k-mers are not found in the dictionary. (For all details, refer to the paper.)

Below a comparison between the dictionary built in Example 2 (not canonical)
and the one just built (Example 3, canonical).

    ./query salmonella_100.index ../data/queries/SRR5833294.10K.fastq.gz
    index size: 10.3981 [MB] (6.36232 [bits/kmer])
    ==== query report:
    num_kmers = 460000
    num_valid_kmers = 459143 (99.8137% of kmers)
    num_positive_kmers = 46 (0.0100187% of valid kmers)
    num_searches = 42/46 (91.3043%)
    num_extensions = 4/46 (8.69565%)
    elapsed = 229.159 millisec / 0.229159 sec / 0.00381932 min / 498.172 ns/kmer

    ./query salmonella_100.canon.index ../data/queries/SRR5833294.10K.fastq.gz
    index size: 11.0657 [MB] (6.77083 [bits/kmer])
    ==== query report:
    num_kmers = 460000
    num_valid_kmers = 459143 (99.8137% of kmers)
    num_positive_kmers = 46 (0.0100187% of valid kmers)
    num_searches = 42/46 (91.3043%)
    num_extensions = 4/46 (8.69565%)
    elapsed = 107.911 millisec / 0.107911 sec / 0.00179852 min / 234.589 ns/kmer

We see that the canonical dictionary is twice as fast as the regular dictionary
for low-hit workloads,
even on this tiny example, for only +0.4 bits/k-mer. -->

Input Files
-----------

LPHash is meant to index k-mers from collections that do not contain duplicates.
These collections can be obtained, for example, by extracting the maximal unitigs of a de Bruijn graph.

To do so, we can use the tool [BCALM2](https://github.com/GATB/bcalm).
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

<!-- #### Datasets
The script `scripts/download_and_preprocess_datasets.sh`
contains all the needed steps to download and pre-process
the datasets that we used in [1]. -->

<!-- Large-Scale Benchmark
---------------------

*Pinus Taeda* ("pine", [GCA_000404065.3](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/065/GCA_000404065.3_Ptaeda2.0/GCA_000404065.3_Ptaeda2.0_genomic.fna.gz)) and *Ambystoma Mexicanum* ("axolotl", [GCA_002915635.2](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/915/635/GCA_002915635.3_AmbMex60DD/GCA_002915635.3_AmbMex60DD_genomic.fna.gz))
are some of the largest genome assemblies, respectively counting
10,508,232,575 and 17,987,935,180 distinct k-mers for k = 31.

After running BCALM2 and UST, we build the indexes as follows.

    ./build ~/DNA_datasets.larger/GCA_000404065.3_Ptaeda2.0_genomic.ust_k31.fa.gz 31 20 -l 6 -c 7 -o pinus.m20.index
    ./build ~/DNA_datasets.larger/GCA_000404065.3_Ptaeda2.0_genomic.ust_k31.fa.gz 31 19 -l 6 -c 7 --canonical-parsing -o pinus.m19.canon.index
    ./build ~/DNA_datasets.larger/GCA_002915635.3_AmbMex60DD_genomic.ust_k31.fa.gz 31 21 -l 6 -c 7 -o axolotl.m21.index
    ./build ~/DNA_datasets.larger/GCA_002915635.3_AmbMex60DD_genomic.ust_k31.fa.gz 31 20 -l 6 -c 7 --canonical-parsing -o axolotl.m20.canon.index

The following table summarizes the space of the dictionaries.

| Dictionary        |Pine       || Axolotl  ||
|:------------------|:---:|:----:|:---:|:---:|
|                   | GB     | bits/k-mer  | GB    | bits/k-mer |
| SSHash, regular   | 13.21  | 10.06       | 22.28 | 9.91       |
| SSHash, canonical | 14.94  | 11.37       | 25.03 | 11.13      |



To query the dictionaries, we use [SRR17023415](https://www.ebi.ac.uk/ena/browser/view/SRR17023415) fastq reads
(23,891,117 reads, each of 150 bases) for the pine,
and [GSM5747680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5747680) multi-line fasta (15,548,160 lines) for the axolotl.

Timings have been collected on an Intel Xeon Platinum 8276L CPU @ 2.20GHz,
using a single thread.

| Dictionary        |Pine       || Axolotl  ||
|:------------------|:---:|:----:|:---:|:---:|
|                   |(>75% hits)||(>86% hits)|
|                   | tot (min) | avg (ns/k-mer) | tot (min) | avg (ns/k-mer) |
| SSHash, regular   | 19.2      | 400            | 4.2       | 269            |
| SSHash, canonical | 14.8      | 310            | 3.2       | 208            |

Below the complete query reports.

    ./query pinus.m20.index ~/DNA_datasets.larger/queries/SRR17023415_1.fastq.gz
    ==== query report:
    num_kmers = 2866934040
    num_valid_kmers = 2866783488 (99.9947% of kmers)
    num_positive_kmers = 2151937575 (75.0645% of valid kmers)
    num_searches = 418897117/2151937575 (19.466%)
    num_extensions = 1733040458/2151937575 (80.534%)
    elapsed = 1146.58 sec / 19.1097 min / 399.933 ns/kmer

    ./query pinus.m19.canon.index ~/DNA_datasets.larger/queries/SRR17023415_1.fastq.gz
    ==== query report:
    num_kmers = 2866934040
    num_valid_kmers = 2866783488 (99.9947% of kmers)
    num_positive_kmers = 2151937575 (75.0645% of valid kmers)
    num_searches = 359426304/2151937575 (16.7025%)
    num_extensions = 1792511271/2151937575 (83.2975%)
    elapsed = 889.779 sec / 14.8297 min / 310.359 ns/kmer

    ./query axolotl.m21.index ~/DNA_datasets.larger/queries/Axolotl.Trinity.CellReports2017.fasta.gz --multiline
    ==== query report:
    num_kmers = 931366757
    num_valid_kmers = 748445346 (80.3599% of kmers)
    num_positive_kmers = 650467884 (86.9092% of valid kmers)
    num_searches = 124008258/650467884 (19.0645%)
    num_extensions = 526459626/650467884 (80.9355%)
    elapsed = 250.173 sec / 4.16955 min / 268.608 ns/kmer

    ./query axolotl.m20.canon.index ~/DNA_datasets.larger/queries/Axolotl.Trinity.CellReports2017.fasta.gz --multiline
    ==== query report:
    num_kmers = 931366757
    num_valid_kmers = 748445346 (80.3599% of kmers)
    num_positive_kmers = 650467884 (86.9092% of valid kmers)
    num_searches = 106220473/650467884 (16.3299%)
    num_extensions = 544247411/650467884 (83.6701%)
    elapsed = 193.871 sec / 3.23119 min / 208.158 ns/kmer -->

Authors
------

[Yoshihiro Shibuya](https://github.com/yhhshb) - <yoshihiro.shibuya@esiee.fr>
[Giulio Ermanno Pibiri](https://jermp.github.io) - <giulioermanno.pibiri@unive.it>
[Antoine Limasset](https://github.com/Malfoy) - <antoine.limasset@univ-lille.fr>

References
-----
* [1] Giulio Ermanno Pibiri and Yoshihiro Shibuya and Antoine Limasset. [Locality-Preserving Minimal Perfect Hashing of k-mers](Submitted).

<!-- ### Datasets

All datasets can be found in `/data2/DNA/lphash_datasets/` on the machine `zipottero.isti.cnr.it`.

Unitigs extracted using bcalm2:

    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 31 -abundance-min 1 -nb-cores 8 -out sal.k31
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 35 -abundance-min 1 -nb-cores 8 -out sal.k35
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 39 -abundance-min 1 -nb-cores 8 -out sal.k39
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 43 -abundance-min 1 -nb-cores 8 -out sal.k43
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 47 -abundance-min 1 -nb-cores 8 -out sal.k47
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 51 -abundance-min 1 -nb-cores 8 -out sal.k51
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 55 -abundance-min 1 -nb-cores 8 -out sal.k55
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 59 -abundance-min 1 -nb-cores 8 -out sal.k59
    ./bcalm -in ~/Salmonella_enterica/Genomes/SAL_CA7616AA.fasta -kmer-size 63 -abundance-min 1 -nb-cores 8 -out sal.k63

    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 31 -abundance-min 1 -nb-cores 8 -out yeast.k31
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 35 -abundance-min 1 -nb-cores 8 -out yeast.k35
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 39 -abundance-min 1 -nb-cores 8 -out yeast.k39
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 43 -abundance-min 1 -nb-cores 8 -out yeast.k43
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 47 -abundance-min 1 -nb-cores 8 -out yeast.k47
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 51 -abundance-min 1 -nb-cores 8 -out yeast.k51
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 55 -abundance-min 1 -nb-cores 8 -out yeast.k55
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 59 -abundance-min 1 -nb-cores 8 -out yeast.k59
    ./bcalm -in ./GCF_000146045.2_R64_genomic.fna.gz -kmer-size 63 -abundance-min 1 -nb-cores 8 -out yeast.k63

    ./bcalm -in ./celegans.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 16 -out celegans.k31
    ./bcalm -in ./celegans.fa.gz -kmer-size 35 -abundance-min 1 -nb-cores 16 -out celegans.k35
    ./bcalm -in ./celegans.fa.gz -kmer-size 39 -abundance-min 1 -nb-cores 16 -out celegans.k39
    ./bcalm -in ./celegans.fa.gz -kmer-size 43 -abundance-min 1 -nb-cores 16 -out celegans.k43
    ./bcalm -in ./celegans.fa.gz -kmer-size 47 -abundance-min 1 -nb-cores 16 -out celegans.k47
    ./bcalm -in ./celegans.fa.gz -kmer-size 51 -abundance-min 1 -nb-cores 16 -out celegans.k51
    ./bcalm -in ./celegans.fa.gz -kmer-size 55 -abundance-min 1 -nb-cores 16 -out celegans.k55
    ./bcalm -in ./celegans.fa.gz -kmer-size 59 -abundance-min 1 -nb-cores 16 -out celegans.k59
    ./bcalm -in ./celegans.fa.gz -kmer-size 63 -abundance-min 1 -nb-cores 16 -out celegans.k63

USTigs computed as follows:

    ./ust -k 31 -i ~/bcalm/build/sal.k31.unitigs.fa
    ./ust -k 35 -i ~/bcalm/build/sal.k35.unitigs.fa
    ./ust -k 39 -i ~/bcalm/build/sal.k39.unitigs.fa
    ./ust -k 43 -i ~/bcalm/build/sal.k43.unitigs.fa
    ./ust -k 47 -i ~/bcalm/build/sal.k47.unitigs.fa
    ./ust -k 51 -i ~/bcalm/build/sal.k51.unitigs.fa
    ./ust -k 55 -i ~/bcalm/build/sal.k55.unitigs.fa
    ./ust -k 59 -i ~/bcalm/build/sal.k59.unitigs.fa
    ./ust -k 63 -i ~/bcalm/build/sal.k63.unitigs.fa
    gzip sal.k31.unitigs.fa.ust.fa
    gzip sal.k35.unitigs.fa.ust.fa
    gzip sal.k39.unitigs.fa.ust.fa
    gzip sal.k43.unitigs.fa.ust.fa
    gzip sal.k47.unitigs.fa.ust.fa
    gzip sal.k51.unitigs.fa.ust.fa
    gzip sal.k55.unitigs.fa.ust.fa
    gzip sal.k59.unitigs.fa.ust.fa
    gzip sal.k63.unitigs.fa.ust.fa

    ./ust -k 31 -i ~/bcalm/build/yeast.k31.unitigs.fa
    ./ust -k 35 -i ~/bcalm/build/yeast.k35.unitigs.fa
    ./ust -k 39 -i ~/bcalm/build/yeast.k39.unitigs.fa
    ./ust -k 43 -i ~/bcalm/build/yeast.k43.unitigs.fa
    ./ust -k 47 -i ~/bcalm/build/yeast.k47.unitigs.fa
    ./ust -k 51 -i ~/bcalm/build/yeast.k51.unitigs.fa
    ./ust -k 55 -i ~/bcalm/build/yeast.k55.unitigs.fa
    ./ust -k 59 -i ~/bcalm/build/yeast.k59.unitigs.fa
    ./ust -k 63 -i ~/bcalm/build/yeast.k63.unitigs.fa
    gzip yeast.k31.unitigs.fa.ust.fa
    gzip yeast.k35.unitigs.fa.ust.fa
    gzip yeast.k39.unitigs.fa.ust.fa
    gzip yeast.k43.unitigs.fa.ust.fa
    gzip yeast.k47.unitigs.fa.ust.fa
    gzip yeast.k51.unitigs.fa.ust.fa
    gzip yeast.k55.unitigs.fa.ust.fa
    gzip yeast.k59.unitigs.fa.ust.fa
    gzip yeast.k63.unitigs.fa.ust.fa

    ./ust -k 31 -i ~/bcalm/build/celegans.k31.unitigs.fa
    ./ust -k 35 -i ~/bcalm/build/celegans.k35.unitigs.fa
    ./ust -k 39 -i ~/bcalm/build/celegans.k39.unitigs.fa
    ./ust -k 43 -i ~/bcalm/build/celegans.k43.unitigs.fa
    ./ust -k 47 -i ~/bcalm/build/celegans.k47.unitigs.fa
    ./ust -k 51 -i ~/bcalm/build/celegans.k51.unitigs.fa
    ./ust -k 55 -i ~/bcalm/build/celegans.k55.unitigs.fa
    ./ust -k 59 -i ~/bcalm/build/celegans.k59.unitigs.fa
    ./ust -k 63 -i ~/bcalm/build/celegans.k63.unitigs.fa
    gzip celegans.k31.unitigs.fa.ust.fa
    gzip celegans.k35.unitigs.fa.ust.fa
    gzip celegans.k39.unitigs.fa.ust.fa
    gzip celegans.k43.unitigs.fa.ust.fa
    gzip celegans.k47.unitigs.fa.ust.fa
    gzip celegans.k51.unitigs.fa.ust.fa
    gzip celegans.k55.unitigs.fa.ust.fa
    gzip celegans.k59.unitigs.fa.ust.fa
    gzip celegans.k63.unitigs.fa.ust.fa -->
