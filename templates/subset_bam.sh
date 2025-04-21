#!/bin/bash

subset-bam_linux --bam '/mnt/c/users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/possorted_genome_bam.bam' --cell-barcodes data/barcodes/Hashtag1_barcodes.tsv --out-bam data/subset_bams/ccrcc_hashtag1.bam --cores 10
 
subset-bam_linux --bam '/mnt/c/users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/possorted_genome_bam.bam' --cell-barcodes data/barcodes/Hashtag2_barcodes.tsv --out-bam data/subset_bams/ccrcc_hashtag2.bam --cores 10

subset-bam_linux --bam '/mnt/c/users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/possorted_genome_bam.bam' --cell-barcodes data/barcodes/Hashtag3_barcodes.tsv --out-bam data/subset_bams/ccrcc_hashtag3.bam --cores 10

subset-bam_linux --bam '/mnt/c/users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/possorted_genome_bam.bam' --cell-barcodes data/barcodes/Hashtag4_barcodes.tsv --out-bam data/subset_bams/ccrcc_hashtag4.bam --cores 10

subset-bam_linux --bam '/mnt/c/users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/possorted_genome_bam.bam' --cell-barcodes data/barcodes/Hashtag5_barcodes.tsv --out-bam data/subset_bams/ccrcc_hashtag5.bam --cores 10

subset-bam_linux --bam '/mnt/c/users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/possorted_genome_bam.bam' --cell-barcodes data/barcodes/Hashtag6_barcodes.tsv --out-bam data/subset_bams/ccrcc_hashtag6.bam --cores 10