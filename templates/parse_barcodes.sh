#!/bin/bash

# --------------------------------------------------
# Shell script to parse and merge cell barcode files
# --------------------------------------------------

# notes:
# - parse cell barcode files to contain unique sample IDs matching BAM files
# - there are 3 samples in HGSOC dataset, and 6 samples in lung dataset

# runtime: seconds

# qsub -V -cwd -l mem_free=10G,h_vmem=20G,h_fsize=100G parse_and_merge_barcodes.sh

# arguments:
# $1: barcodes file
# $2: sample ID

# -----------------------------------
# start runtime
#start=`date +%s`
# -----------------------------------


# mkdir -p barcodes_merged

#cp data/barcodes/$1_barcodes.tsv data/barcodes_merged/barcodes_$7.tsv

# add unique sample IDs to cell barcodes for each sample
# sed -i "s|\([A-Z]\+\)\-1|\1\-$7|g" data/barcodes_merged/barcodes_$7.tsv
sed "s|\([A-Z]\+\)\-1|\1\-K${2}|g" $1 > ${1}_rn.tsv


# -----------------------------------
# end runtime
# end=`date +%s`
# runtime=`expr $end - $start`

# save runtime
# mkdir -p $1/parse_and_merge_barcodes
# echo runtime: $runtime seconds > $1/parse_and_merge_barcodes/runtime_parse_and_merge_barcodes.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
# mkdir -p $2/parse_and_merge_barcodes
# date > $2/parse_and_merge_barcodes/timestamp_parse_and_merge_barcodes.txt
# -----------------------------------

