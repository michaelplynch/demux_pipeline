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
# $1: tsv 1
# $2: tsv 2
# $3: tsv 3
# $4: tsv 4
# $5: tsv 5
# $6: tsv 6
# $7: short sample ID 1
# $8: short sample ID 2
# $9: short sample ID 3
# $10: short sample ID 4
# $11: short sample ID 5
# $12: short sample ID 6

# -----------------------------------
# start runtime
#start=`date +%s`
# -----------------------------------


mkdir -p data/barcodes_merged

cp data/barcodes/$1_barcodes.tsv data/barcodes_merged/barcodes_$7.tsv
cp data/barcodes/$2_barcodes.tsv data/barcodes_merged/barcodes_$8.tsv
cp data/barcodes/$3_barcodes.tsv data/barcodes_merged/barcodes_$9.tsv
cp data/barcodes/$4_barcodes.tsv data/barcodes_merged/barcodes_${10}.tsv
cp data/barcodes/$5_barcodes.tsv data/barcodes_merged/barcodes_${11}.tsv
cp data/barcodes/$6_barcodes.tsv data/barcodes_merged/barcodes_${12}.tsv

# add unique sample IDs to cell barcodes for each sample
sed -i "s|\([A-Z]\+\)\-1|\1\-$7|g" data/barcodes_merged/barcodes_$7.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-$8|g" data/barcodes_merged/barcodes_$8.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-$9|g" data/barcodes_merged/barcodes_$9.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${10}|g" data/barcodes_merged/barcodes_${10}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${11}|g" data/barcodes_merged/barcodes_${11}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${12}|g" data/barcodes_merged/barcodes_${12}.tsv

# merge files
cat data/barcodes_merged/barcodes_$7.tsv data/barcodes_merged/barcodes_$8.tsv data/barcodes_merged/barcodes_$9.tsv data/barcodes_merged/barcodes_${10}.tsv data/barcodes_merged/barcodes_${11}.tsv data/barcodes_merged/barcodes_${12}.tsv > \
data/barcodes_merged/barcodes_merged.tsv


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

