THIS_PATH=$(echo $PWD)
LPHASH_DIR=$(dirname $THIS_PATH)
UNITIGS_FOLDER=$LPHASH_DIR/data/unitigs_stitched
QUERY_FOLDER=$LPHASH_DIR/data/queries
MPHF_FOLDER=$LPHASH_DIR/mphfs
RESULTS_FOLDER=$LPHASH_DIR/results

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

## -----------------------------------------------------------------------------

THREADS=1
C=4.0

M="15"
QUERY="$QUERY_FOLDER/salmonella_enterica.fasta.gz"
# for K in 31 35 39 43 47 51 55 59 63 ; do
for K in 63 ; do
    UNITIGS=$UNITIGS_FOLDER/"se.ust.k$K.fa.gz"
    LPMPHF="$MPHF_FOLDER/se.k$K.m$M.lphash.bin"
    PTMPHF="$MPHF_FOLDER/se.ust.k$K.pthash.bin"
    $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
    $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
    $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c $C >> $PTBUILD_RESULTS
    $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
done

# M=16
# QUERY="$QUERY_FOLDER/saccaromyces_cerevisae.fasta.gz"
# for K in 31 35 39 43 47 51 55 59 63 ; do
#     local UNITIGS=$UNITIGS_FOLDER/"sc.ust.k$K.fa.gz"
#     local LPMPHF="$MPHF_FOLDER/sc.k$K.m$M.lphash.bin"
#     local PTMPHF="$MPHF_FOLDER/sc.ust.k$K.pthash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
# done

# M=17
# QUERY="$QUERY_FOLDER/celegans.fasta.gz"
# for K in 31 35 39 43 47 51 55 59 63 ; do
#     local UNITIGS=$UNITIGS_FOLDER/"celegans.ust.k$K.fa.gz"
#     local LPMPHF="$MPHF_FOLDER/celegans.k$K.m$M.lphash.bin"
#     local PTMPHF="$MPHF_FOLDER/celegans.ust.k$K.pthash.bin"
#     $LPBUILD $UNITIGS $K $M -t $THREADS -o $LPMPHF -d "tmp/" >> $LPBUILD_RESULTS
#     $LPQUERY $LPMPHF $QUERY >> $LPQUERY_RESULTS
#     $PTBUILD $UNITIGS $K -t $THREADS -p $PTMPHF -d "tmp/" -c 4.0 >> $PTBUILD_RESULTS
#     $PTQUERY $QUERY $K -p $PTMPHF >> $PTQUERY_RESULTS
# done

