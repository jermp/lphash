#!/bin/sh

THIS_PATH=$(echo $PWD)
DATASETS_PATH="/data2/DNA/lphash_datasets"
BINARIES_PATH="/data2/DNA/lphash_binaries"
LPHASH_DIR="/home/shibuya/lphash"
BUILD_DIR=$LPHASH_DIR/build
COMPILE_OPTIONS=$LPHASH_DIR/"include/compile_constants.tpd"
UNITIGS_FOLDER=$DATASETS_PATH
QUERY_FOLDER=$LPHASH_DIR/data/queries
MPHF_FOLDER=$LPHASH_DIR/mphfs
MPHF_FOLDER=$BINARIES_PATH
RESULTS_FOLDER=$LPHASH_DIR/results
TMP_FOLDER="/data2/DNA/tmp_dir"

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

THREADS="4"
C="5.0"
QUERY="/data2/DNA/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz"

## -----------------------------------------------------------------------------

echo "// typedef kmer128_t kmer_t;\ntypedef uint64_t kmer_t;" > $COMPILE_OPTIONS
cd $BUILD_DIR
make -j
cd $THIS_PATH

## -----------------------------------------------------------------------------

for K in 31 ; do
    for M in 15 ; do
        UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 31 ; do
    for M in 16 ; do
        UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 31 ; do
    for M in 20 ; do
        UNITIGS=$UNITIGS_FOLDER/"cod.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/cod.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 31 ; do
    for M in 20 ; do
        UNITIGS=$UNITIGS_FOLDER/"kestrel.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/kestrel.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 31 ; do
    for M in 21 ; do
        UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

## -----------------------------------------------------------------------------

echo "typedef kmer128_t kmer_t;\n// typedef uint64_t kmer_t;" > $COMPILE_OPTIONS
cd $BUILD_DIR
make -j
cd $THIS_PATH

## Yeast -----------------------------------------------------------------------------

for K in 35 ; do
    for M in 15 ; do
        UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 39 43 47 51 ; do
    for M in 16 ; do
        UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 55 59 63 ; do
    for M in 18 ; do
        UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/yeast.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

## Celegans -------------------------------------------------------------------------------------------------

for K in 35 39 ; do
    for M in 18 ; do
        UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 43 47 51 55 59 63 ; do
    for M in 20 ; do
        UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

## Cod --------------------------------------------------------------------------------------------------

for K in 35 ; do
    for M in 20 ; do
        UNITIGS=$UNITIGS_FOLDER/"cod.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/cod.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 39 43 47 ; do
    for M in 22 ; do
        UNITIGS=$UNITIGS_FOLDER/"cod.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/cod.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 51 55 59 63 ; do
    for M in 24 ; do
        UNITIGS=$UNITIGS_FOLDER/"cod.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/cod.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

## Kestrel ---------------------------------------------------------------------------------------------------

for K in 35 ; do
    for M in 20 ; do
        UNITIGS=$UNITIGS_FOLDER/"kestrel.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/kestrel.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 39 43 47 ; do
    for M in 22 ; do
        UNITIGS=$UNITIGS_FOLDER/"kestrel.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/kestrel.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 51 55 59 63 ; do
    for M in 24 ; do
        UNITIGS=$UNITIGS_FOLDER/"kestrel.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/kestrel.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

## Human -------------------------------------------------------------------------------------------------------

for K in 35 ; do
    for M in 21 ; do
        UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 39 43 ; do
    for M in 23 ; do
        UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 47 51 ; do
    for M in 26 ; do
        UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

for K in 55 59 63 ; do
    for M in 28 ; do
        UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
        LPMPHF="$MPHF_FOLDER/human.k$K.m$M.lphash.bin"
        $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS
        $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
        $LPBUILD_ALT $UNITIGS $K $M -t $THREADS -o $LPMPHF -d $TMP_FOLDER -c $C >> $LPBUILD_RESULTS_ALT
        $LPQUERY_ALT $LPMPHF $QUERY >> $LPQUERY_RESULTS_ALT
    done
done

## PTHash and BBHash ----------------------------------------------------------------------------------------------------

K=63
UNITIGS=$UNITIGS_FOLDER/"yeast.k$K.unitigs.fa.ust.fa.gz"
PTMPHF="$MPHF_FOLDER/yeast.k$K.pthash.bin"
BBMPHF="$MPHF_FOLDER/yeast.k$K.bbhash.bin"
$PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d $TMP_FOLDER -c $C >> $PTBUILD_RESULTS
$PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS

K=63
UNITIGS=$UNITIGS_FOLDER/"celegans.k$K.unitigs.fa.ust.fa.gz"
PTMPHF="$MPHF_FOLDER/celegans.k$K.pthash.bin"
BBMPHF="$MPHF_FOLDER/celegans.k$K.bbhash.bin"
$PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d $TMP_FOLDER -c $C >> $PTBUILD_RESULTS
$PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS

K=63
UNITIGS=$UNITIGS_FOLDER/"cod.k$K.unitigs.fa.ust.fa.gz"
PTMPHF="$MPHF_FOLDER/cod.k$K.pthash.bin"
BBMPHF="$MPHF_FOLDER/cod.k$K.bbhash.bin"
$PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d $TMP_FOLDER -c $C >> $PTBUILD_RESULTS
$PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS

K=63
UNITIGS=$UNITIGS_FOLDER/"kestrel.k$K.unitigs.fa.ust.fa.gz"
PTMPHF="$MPHF_FOLDER/kestrel.k$K.pthash.bin"
BBMPHF="$MPHF_FOLDER/kestrel.k$K.bbhash.bin"
$PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d $TMP_FOLDER -c $C >> $PTBUILD_RESULTS
$PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS

K=63
UNITIGS=$UNITIGS_FOLDER/"human.k$K.unitigs.fa.ust.fa.gz"
PTMPHF="$MPHF_FOLDER/human.k$K.pthash.bin"
BBMPHF="$MPHF_FOLDER/human.k$K.bbhash.bin"
$PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -b $BBMPHF -d $TMP_FOLDER -c $C >> $PTBUILD_RESULTS
$PTQUERY $QUERY $K -p $PTMPHF -b $BBMPHF >> $PTQUERY_RESULTS