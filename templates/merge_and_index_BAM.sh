#!/bin/bash

# ------------------------------------------------
# Shell script to merge and index parsed BAM files
# ------------------------------------------------

# notes:
# - BAM files from Cell Ranger are already position sorted, so do not need to sort
# - there are 3 samples in HGSOC dataset, and 6 samples in lung dataset

# runtime: ~4-6 hours

# qsub -V -cwd -l mem_free=10G,h_vmem=20G,h_fsize=100G merge_and_index_BAM.sh

# arguments:
# $1: input bams
# $2: merged bam
# $3  outdir

# -----------------------------------
# start runtime
# start=`date +%s`
# -----------------------------------


mkdir -p $3

# merge BAM files
samtools merge -@ 4 $3/$2 \
$1

# index merged BAM
samtools index $3/$2


# -----------------------------------
# end runtime
# end=`date +%s`
# runtime=`expr $end - $start`

# save runtime
# mkdir -p $1/merge_and_index_BAM
# echo runtime: $runtime seconds > $1/merge_and_index_BAM/runtime_merge_and_index_BAM.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
# mkdir -p $2/merge_and_index_BAM
# date > $2/merge_and_index_BAM/timestamp_merge_and_index_BAM.txt
# -----------------------------------

