# $1 barcodes file
# $2 key
# $3 number of cells for this hashtag

args <- commandArgs(trailingOnly = TRUE)
barcodes<-args[1]
key<-args[2]
n<-args[3]
hashtag<-args[4]

library(tidyverse)
library(magrittr)
#barcodes<-'/data/projects/demuxSNP/data/processed/10x_cellline_2x/raji/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji_count_sample_barcodes.csv'
print(barcodes)

bcs<-read_table(barcodes, col_names = "barcode")

print(length(bcs$barcode))
print(key);print(hashtag);print(n);print(barcodes)
set.seed(1)
bcs_sub<-bcs$barcode[sample(seq_along(bcs$barcode),n)]
write_tsv(tibble(bcs_sub),file=paste("Hashtag",hashtag,"_", key,"_",n,"_sub.tsv",sep=""),col_names=FALSE)