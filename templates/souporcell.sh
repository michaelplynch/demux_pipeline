#!/bin/bash

## record code used for souporcell command

sudo singularity exec -B /home/mlynch/souporcell /home/mlynch/mybin/souporcell_latest.sif souporcell_pipeline.py -i doublet_bam/bam_merged_doublets_ccrcc_40pc.bam.bam -b barcodes_merged_ccrcc_40pc.tsv -f genome.fa -t 10 -o souporcell_results_40pc -k 6 --skip_remap SKIP_REMAP --common_variants genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf

time