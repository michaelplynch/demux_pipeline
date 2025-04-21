#!/bin/bash

# arguments:
# $1: input bam
# $2: output bam
# $3: short sample ID

./scripts/parse_BAM_files.sh ccrcc_hashtag1.bam ccrcc_hashtag1_parsed.bam K1

./scripts/parse_BAM_files.sh ccrcc_hashtag2.bam ccrcc_hashtag2_parsed.bam K2

./scripts/parse_BAM_files.sh ccrcc_hashtag3.bam ccrcc_hashtag3_parsed.bam K3

./scripts/parse_BAM_files.sh ccrcc_hashtag4.bam ccrcc_hashtag4_parsed.bam K4

./scripts/parse_BAM_files.sh ccrcc_hashtag5.bam ccrcc_hashtag5_parsed.bam K5

./scripts/parse_BAM_files.sh ccrcc_hashtag6.bam ccrcc_hashtag6_parsed.bam K6