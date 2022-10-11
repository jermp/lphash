#!/bin/sh

THIS_PATH=$(echo $PWD)
DATASETS_PATH="/data2/DNA/lphash_datasets"
BINARIES_PATH="/data2/DNA/lphash_binaries"
LPHASH_DIR="/home/shibuya/lphash"
#$(dirname $THIS_PATH)
BUILD_DIR=$LPHASH_DIR/build
COMPILE_OPTIONS=$LPHASH_DIR/"include/compile_constants.tpd"
UNITIGS_FOLDER=$DATASETS_PATH #$LPHASH_DIR/data/unitigs_stitched
QUERY_FOLDER=$LPHASH_DIR/data/queries
# MPHF_FOLDER=$LPHASH_DIR/mphfs
MPHF_FOLDER=$BINARIES_PATH
RESULTS_FOLDER=$LPHASH_DIR/results
TMP_FOLDER="/data2/DNA/tmp_dir"
# $LPHASH_DIR/tmp

LPBUILD=$LPHASH_DIR/build/build
LPQUERY=$LPHASH_DIR/build/query
LPBUILD_ALT=$LPHASH_DIR/build/build_alt
LPQUERY_ALT=$LPHASH_DIR/build/query_alt
PTBUILD=$LPHASH_DIR/build/ptbb_build
PTQUERY=$LPHASH_DIR/build/ptbb_query

LPBUILD_RESULTS=$LPHASH_DIR/results/lphash_build.csv
LPQUERY_RESULTS=$LPHASH_DIR/results/lphash_query.csv
LPBUILD_RESULTS_ALT=$LPHASH_DIR/results/lphash_build_alt.csv
LPQUERY_RESULTS_ALT=$LPHASH_DIR/results/lphash_query_alt.csv
PTBUILD_RESULTS=$LPHASH_DIR/results/pthash_build.csv
PTQUERY_RESULTS=$LPHASH_DIR/results/pthash_query.csv

## -----------------------------------------------------------------------------

mkdir -p $MPHF_FOLDER;
mkdir -p $RESULTS_FOLDER;
mkdir -p $TMP_FOLDER;

## -----------------------------------------------------------------------------

if [ ! -f $LPBUILD_RESULTS ]
then
    echo "unitigs file,k,m,fraction of non-unique minimizers,estimated epsilon,true epsilon,alpha,space" > $LPBUILD_RESULTS
fi

if [ ! -f $LPQUERY_RESULTS ]
then
    echo "input file,lphash file,total k-mers,barebone streaming time,barebone random time,full streaming time,full random time" > $LPQUERY_RESULTS
fi

if [ ! -f $LPBUILD_RESULTS_ALT ]
then
    echo "unitigs file,k,m,fraction of non-unique minimizers,estimated epsilon,true epsilon,alpha,space" > $LPBUILD_RESULTS_ALT
fi

if [ ! -f $LPQUERY_RESULTS_ALT ]
then
    echo "input file,lphash (alt) file,total k-mers,barebone streaming time,barebone random time,full streaming time,full random time" > $LPQUERY_RESULTS_ALT
fi

if [ ! -f $PTBUILD_RESULTS ]
then
    echo "unitigs file,k,total k-mers,pthash size,pthash space,bbhash size,bbhash space" > $PTBUILD_RESULTS
fi

if [ ! -f $PTQUERY_RESULTS ]
then
    echo "input file,k,pthash file,pt time,bbhash file,bb time" > $PTQUERY_RESULTS
fi

## -----------------------------------------------------------------------------

THREADS=4
C=5.0

## -----------------------------------------------------------------------------

# echo "// typedef kmer128_t kmer_t;\ntypedef uint64_t kmer_t;" > $COMPILE_OPTIONS
# cd $BUILD_DIR
# make -j
# cd $THIS_PATH

## -----------------------------------------------------------------------------

# M=15
# QUERY="$QUERY_FOLDER/salmonella_enterica.fasta.gz"
# for K in 31 ; do
#     UNITIGS=$UNITIGS_FOLDER/"sal.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/sal.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/sal.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/sal.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

# M=16
# # QUERY="$QUERY_FOLDER/saccaromyces_cerevisae.fasta.gz"
# for K in 31 ; do
#     UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/yeast.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/yeast.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

# M=20
# # QUERY="$QUERY_FOLDER/celegans.fasta.gz"
# for K in 31 ; do
#     UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/celegans.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/celegans.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

# M=20
# # QUERY="$QUERY_FOLDER/cod.fasta.gz"
# for K in 31 ; do
#     UNITIGS=$UNITIGS_FOLDER/"cod.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/cod.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/cod.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/cod.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

# M=27
# # QUERY="$QUERY_FOLDER/human.fasta.gz"
# for K in 31 ; do
#     UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/human.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/human.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

## -----------------------------------------------------------------------------

echo "typedef kmer128_t kmer_t;\n// typedef uint64_t kmer_t;" > $COMPILE_OPTIONS
cd $BUILD_DIR
make -j
cd $THIS_PATH

## -----------------------------------------------------------------------------

# M=15
# QUERY="$QUERY_FOLDER/salmonella_enterica.fasta.gz"
# for K in 35 39 43 47 51 55 59 63 ; do
#     UNITIGS=$UNITIGS_FOLDER/"sal.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/sal.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/sal.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/sal.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

# M=16
# # QUERY="$QUERY_FOLDER/saccaromyces_cerevisae.fasta.gz"
# for K in 35 39 43 47 51 55 59 63 ; do
#     UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/yeast.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/yeast.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

# M=20
# # QUERY="$QUERY_FOLDER/celegans.fasta.gz"
# for K in 35 39 43 47 51 55 59 63 ; do
#     UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
#     PTMPHF="$MPHF_FOLDER/celegans.k$K.pthash.bin"
#     BBMPHF="$MPHF_FOLDER/celegans.k$K.bbhash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS
# done

# M=24
# QUERY="$QUERY_FOLDER/kestrel.fasta.gz"
# QUERY="$QUERY_FOLDER/salmonella_enterica.fasta.gz"
# for K in 43 47 51 55 59 63 ; do
#     UNITIGS=$UNITIGS_FOLDER/"kestrel.k$K.unitigs.fa.ust.fa.gz"
#     LPMPHF="$MPHF_FOLDER/kestrel.k$K.m$M.lphash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#     $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
# done

# K=63
# PTMPHF="$MPHF_FOLDER/kestrel.k$K.pthash.bin"
# BBMPHF="$MPHF_FOLDER/kestrel.k$K.bbhash.bin"
# $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
# $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS

# M=20
# QUERY="$QUERY_FOLDER/cod.fasta.gz"
# QUERY="$QUERY_FOLDER/salmonella_enterica.fasta.gz"
# for K in 43 47 51 55 59 63 ; do
#     for M in 20 21 22 23 24 25 ; do
#         UNITIGS=$UNITIGS_FOLDER/"cod.k$K.unitigs.fa.ust.fa.gz"
#         LPMPHF="$MPHF_FOLDER/cod.k$K.m$M.lphash.bin"
#         $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#         $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#         $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
#         $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
#     done
# done

# K=63
# PTMPHF="$MPHF_FOLDER/cod.k$K.pthash.bin"
# BBMPHF="$MPHF_FOLDER/cod.k$K.bbhash.bin"
# $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
# $PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS

M=27
QUERY="/data2/DNA/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz"
for K in 43 47 51 55 59 63 ; do
    UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
    LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS_ALT
    $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
done

K=63
PTMPHF="$MPHF_FOLDER/human.k$K.pthash.bin"
BBMPHF="$MPHF_FOLDER/human.k$K.bbhash.bin"
$PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
$PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS

## -----------------------------------------------------------------------------


