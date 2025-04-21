#!/bin/bash

# ------------------------------------------------------------------
# Shell script to parse BAM files to add sample IDs to cell barcodes
# ------------------------------------------------------------------

# notes:
# - convert BAM to SAM, parse SAM, then convert back to BAM to save space
# - syntax to search and replace: sed -i "s/regexp/replacement/g"
# - using double quotes allows variables inside sed expression
# - using alternative separator | allows slashes in variable inside sed expression
# - regular expression matches cell barcode with format "CB:Z:AACTTTCAGCGCTCCA-1"
# - we replace the sample suffix "-1" with a unique sample ID, e.g. "-X2"

# runtime: ~1-2 hours

# qsub -V -cwd -l mem_free=10G,h_vmem=20G,h_fsize=100G parse_BAM_files.sh

# arguments:
# $1: input bam
# $2: output bam
# $3: short sample ID
# $4: outdir

# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------

mkdir -p $4

samtools view -h $1 | \
sed "s|\(CB\:Z\:[A-Z]\+\)\-1|\1\-$3|g" | \
samtools view -bo $4/$2

# -----------------------------------
# end runtime
# end=`date +%s`
# runtime=`expr $end - $start`

# save runtime
# mkdir -p parse_BAM_files
# echo runtime: $runtime seconds > parse_BAM_files/runtime_parse_BAM_files_$3.txt
# -----------------------------------