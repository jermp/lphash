#!/bin/bash

function debian_realpath() {
    f=$@;
    if [ -d "$f" ]; then
        base="";
        dir="$f";
    else
        base="/$(basename "$f")";
        dir=$(dirname "$f");
    fi
    dir=$(cd "$dir" && /bin/pwd);
    echo "$dir$base";
}

script_name=$0
script_full_path=$(dirname $(debian_realpath "$0"))
lphash_full_path=$(dirname "$script_full_path")

unitigs_folder="$lphash_full_path/data/unitigs_stitched"
queries_folder="$lphash_full_path/data/queries"

mkdir -p $unitigs_folder
mkdir -p $queries_folder 

## queries

ECOLI1=ecoli1.fasta.gz
SRR5833294=SRR5833294.10K.fastq.gz
SALMONELLA=salmonella_enterica.fasta.gz

wget "https://github.com/yhhshb/test_datasets/raw/main/escherichia_coli/genome/$ECOLI1"
wget "https://github.com/yhhshb/test_datasets/raw/main/human/reads/$SRR5833294"
wget "https://github.com/yhhshb/test_datasets/blob/main/salmonella_enterica/genome/$SALMONELLA"

mv $ECOLI1 $queries_folder/$ECOLI1
mv $SRR5833294 $queries_folder/$SRR5833294
mv $SALMONELLA $queries_folder/$SALMONELLA

## input k-mer sets

SAK31=se.ust.k31.fa.gz
SAK47=se.ust.k47.fa.gz
SAK63=se.ust.k63.fa.gz

wget "https://github.com/yhhshb/test_datasets/raw/main/salmonella_enterica/unitigs/$SAK31"
wget "https://github.com/yhhshb/test_datasets/raw/main/salmonella_enterica/unitigs/$SAK47"
wget "https://github.com/yhhshb/test_datasets/raw/main/salmonella_enterica/unitigs/$SAK63"

mv $SAK31 $unitigs_folder/$SAK31
mv $SAK47 $unitigs_folder/$SAK47
mv $SAK63 $unitigs_folder/$SAK63
