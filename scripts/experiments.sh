#!/bin/sh

THIS_PATH=$(echo $PWD)
DATASETS_PATH="/data2/DNA/lphash_datasets"
LPHASH_DIR="/home/shibuya/lphash"
#$(dirname $THIS_PATH)
BUILD_DIR=$LPHASH_DIR/build
COMPILE_OPTIONS=$LPHASH_DIR/"include/compile_constants.tpd"
UNITIGS_FOLDER=$DATASETS_PATH #$LPHASH_DIR/data/unitigs_stitched
QUERY_FOLDER=$LPHASH_DIR/data/queries
MPHF_FOLDER=$LPHASH_DIR/mphfs
RESULTS_FOLDER=$LPHASH_DIR/results
TMP_FOLDER=$LPHASH_DIR/tmp

LPBUILD=$LPHASH_DIR/build/build
LPQUERY=$LPHASH_DIR/build/query
PTBUILD=$LPHASH_DIR/build/ptbb_build
PTQUERY=$LPHASH_DIR/build/ptbb_query

LPBUILD_RESULTS=$LPHASH_DIR/results/lphash_build.csv
LPQUERY_RESULTS=$LPHASH_DIR/results/lphash_query.csv
PTBUILD_RESULTS=$LPHASH_DIR/results/pthash_build.csv
PTQUERY_RESULTS=$LPHASH_DIR/results/pthash_query.csv

## -----------------------------------------------------------------------------

mkdir -p $MPHF_FOLDER;
mkdir -p $RESULTS_FOLDER;
mkdir -p $TMP_FOLDER;

THREADS=4
C=5.0

## -----------------------------------------------------------------------------

echo "// typedef kmer128_t kmer_t;\ntypedef uint64_t kmer_t;" > $COMPILE_OPTIONS
cd BUILD_DIR
make -j
cd THIS_PATH

## -----------------------------------------------------------------------------

M=15
QUERY="$QUERY_FOLDER/salmonella_enterica.fasta.gz"
for K in 31 ; do
    UNITIGS=$UNITIGS_FOLDER/"sal.k$K.unitigs.fa.ust.fa.gz"
    LPMPHF="$MPHF_FOLDER/sal.k$K.m$M.lphash.bin"
    PTMPHF="$MPHF_FOLDER/sal.k$K.pthash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
    $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
done

M=16
# QUERY="$QUERY_FOLDER/saccaromyces_cerevisae.fasta.gz"
for K in 31 ; do
    UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
    LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
    PTMPHF="$MPHF_FOLDER/yeast.k$K.pthash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
    $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
done

M=17
# QUERY="$QUERY_FOLDER/celegans.fasta.gz"
for K in 31 ; do
    UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
    LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
    PTMPHF="$MPHF_FOLDER/celegans.k$K.pthash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
    $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
done

# M=27
# # QUERY="$QUERY_FOLDER/human.fasta.gz"
# for K in 31 ; do
#     UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/human.k$K.pthash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
# done

## -----------------------------------------------------------------------------

echo "typedef kmer128_t kmer_t;\n// typedef uint64_t kmer_t;" > $COMPILE_OPTIONS
cd BUILD_DIR
make -j
cd THIS_PATH

## -----------------------------------------------------------------------------

M=15
QUERY="$QUERY_FOLDER/salmonella_enterica.fasta.gz"
for K in 35 39 43 47 51 55 59 63 ; do
    UNITIGS=$UNITIGS_FOLDER/"sal.k$K.unitigs.fa.ust.fa.gz"
    LPMPHF="$MPHF_FOLDER/sal.k$K.m$M.lphash.bin"
    PTMPHF="$MPHF_FOLDER/sal.k$K.pthash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
    $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
done

M=16
# QUERY="$QUERY_FOLDER/saccaromyces_cerevisae.fasta.gz"
for K in 35 39 43 47 51 55 59 63 ; do
    UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
    LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
    PTMPHF="$MPHF_FOLDER/yeast.k$K.pthash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
    $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
done

M=17
# QUERY="$QUERY_FOLDER/celegans.fasta.gz"
for K in 35 39 43 47 51 55 59 63 ; do
    UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
    LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
    PTMPHF="$MPHF_FOLDER/celegans.k$K.pthash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
    $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
done

# M=27
# # QUERY="$QUERY_FOLDER/celegans.fasta.gz"
# for K in 35 39 43 47 51 55 59 63 ; do
#     UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/human.k$K.pthash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
# done

## -----------------------------------------------------------------------------


