ZENODO_TOKEN=$1

THIS_PATH=$(echo $PWD)
DATASETS_PATH="/data2/DNA/lphash_datasets"
UNITIGS_FOLDER=$DATASETS_PATH

for K in 31 35 39 43 47 51 55 59 63 ; do
    FNAME="sal.k$K.unitigs.fa.ust.fa.gz"
    UNITIGS=$UNITIGS_FOLDER/$FNAME
    curl --upload-file ${UNITIGS} "https://zenodo.org/api/files/327936a5-d053-4dd4-897f-39bfaa6fad51/${FNAME}?access_token=${ZENODO_TOKEN}"
done

for K in 31 35 39 43 47 51 55 59 63 ; do
    FNAME="yeast.k$K.unitigs.fa.ust.fa.gz"
    UNITIGS=$UNITIGS_FOLDER/$FNAME
    curl --upload-file ${UNITIGS} "https://zenodo.org/api/files/327936a5-d053-4dd4-897f-39bfaa6fad51/${FNAME}?access_token=${ZENODO_TOKEN}"
done

for K in 31 35 39 43 47 51 55 59 63 ; do
    FNAME="celegans.k$K.unitigs.fa.ust.fa.gz"
    UNITIGS=$UNITIGS_FOLDER/$FNAME
    curl --upload-file ${UNITIGS} "https://zenodo.org/api/files/327936a5-d053-4dd4-897f-39bfaa6fad51/${FNAME}?access_token=${ZENODO_TOKEN}"
done

for K in 31 35 39 43 47 51 55 59 63 ; do
    FNAME="cod.k$K.unitigs.fa.ust.fa.gz"
    UNITIGS=$UNITIGS_FOLDER/$FNAME
    curl --upload-file ${UNITIGS} "https://zenodo.org/api/files/327936a5-d053-4dd4-897f-39bfaa6fad51/${FNAME}?access_token=${ZENODO_TOKEN}"
done

for K in 31 35 39 43 47 51 55 59 63 ; do
    FNAME="kestrel.k$K.unitigs.fa.ust.fa.gz"
    UNITIGS=$UNITIGS_FOLDER/$FNAME
    curl --upload-file ${UNITIGS} "https://zenodo.org/api/files/327936a5-d053-4dd4-897f-39bfaa6fad51/${FNAME}?access_token=${ZENODO_TOKEN}"
done

for K in 31 35 39 43 47 51 55 59 63 ; do
    FNAME="human.k$K.unitigs.fa.ust.fa.gz"
    UNITIGS=$UNITIGS_FOLDER/$FNAME
    curl --upload-file ${UNITIGS} "https://zenodo.org/api/files/327936a5-d053-4dd4-897f-39bfaa6fad51/${FNAME}?access_token=${ZENODO_TOKEN}"
done