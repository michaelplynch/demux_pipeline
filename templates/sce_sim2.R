library(demuxSNPpaperfigures)
library(dittoSeq)
library(ggpubr)
library(gridExtra)
library(tictoc)
library(SingleCellExperiment)
library(Matrix)
library(tidyverse)
library(demuxSNP)


args <- commandArgs(trailingOnly = TRUE)
tenx_path<-args[1]#c('/data/projects/yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/filtered_feature_bc_matrix')
barcodes_path<-args[2]#c("/data/projects/demuxSNP/nextflow/ccrcc_out/barcodes_merged_e.tsv")
lookup_path<-args[3]#c("/data/projects/demuxSNP/nextflow/ccrcc_out/lookup_table_doublets_pbmc_e_25pc.tsv")
key<-args[4]#c("e")

matrix<-readMM(paste(tenx_path,"/matrix.mtx.gz",sep=""))
barcodes<-read.table(paste(tenx_path,"/barcodes.tsv.gz",sep=""))
features<-read.table(paste(tenx_path, "/features.tsv.gz",sep=""))

rownames(matrix)<-features$V2
colnames(matrix)<-barcodes$V1

rna<-matrix[features$V3=="Gene",]
hto<-matrix[grep("Hashtag",features$V1),]

sub_barcodes<-read.table(barcodes_path)
lookup<-read_tsv(lookup_path)

rna<-rna[,substr(colnames(rna),start=1,stop = 16) %in% substr(sub_barcodes$V1, start=1,stop=16)]
hto<-hto[,substr(colnames(hto),start=1,stop = 16) %in% substr(sub_barcodes$V1, start=1,stop=16)]

new_bc<-sub_barcodes$V1[match(substr(colnames(rna),1,16),substr(sub_barcodes$V1,1,16))]
colnames(rna)<-new_bc
colnames(hto)<-new_bc
rna1<-rna[,lookup$original]
rna2<-rna[,lookup$replacement]
print(dim(rna1));print(dim(rna2))
rna_doub<-rna2+rna1
rna_sing<-rna[,!c(new_bc %in% union(lookup$original,lookup$replacement))]
rna_full<-cbind(rna_sing,rna_doub)

# % rna per doublet cell
pc_rna<-colSums(rna1)/(colSums(rna1)+colSums(rna2))

hto1<-hto[,lookup$original]
hto2<-hto[,lookup$replacement]
print(dim(hto1));print(dim(hto2))
hto_doub<-hto2+hto1
hto_sing<-hto[,!c(new_bc %in% union(lookup$original,lookup$replacement))]
hto_full<-cbind(hto_sing,hto_doub)

sce<-SingleCellExperiment(list(counts=rna_full))
hto_sce<-SingleCellExperiment(list(counts=hto_full))
altExp(sce,"HTO")<-hto_sce
mainExpName(sce)<-"RNA"
head(colnames(sce))

sce$new_bc<-sub_barcodes$V1[match(substr(colnames(sce),1,16),substr(sub_barcodes$V1,1,16))]

sce$truth<-substr(sce$new_bc,start=18,stop=19)
sce$truth[sce$new_bc %in% lookup$replacement[substr(lookup$original,start=18,stop=19)!=substr(lookup$replacement,start=18,stop=19)]]<-"Doublet"
table(sce$truth)

sce$truthalldoub<-substr(sce$new_bc,start=18,stop=19)
sce$truthalldoub[sce$new_bc %in% lookup$replacement]<-"Doublet"
table(sce$truthalldoub)

sce$truthfull<-substr(sce$new_bc,start=18,stop=19)
order<-na.omit(match(substr(sce$new_bc,1,16),substr(lookup$replacement,1,16)))
sce$truthfull[colnames(sce) %in% lookup$replacement]<-paste(substr(lookup$original,18,19),substr(lookup$replacement,18,19),sep="")[order]
table(sce$truthfull)

rownames(hto)
sig<-matrix(0,dim(hto)[1],dim(hto)[2])
colnames(sig)<-colnames(hto)
t<-substr(new_bc,18,19)
truth<-gsub("K","Hashtag",t)
for (i in seq_along(rownames(hto))) {
  hashtag<-rownames(hto)[i]
  sig[i,truth==hashtag]<-1
}

sig1<-sig[,lookup$original]
sig2<-sig[,lookup$replacement]
print(dim(sig1));print(dim(sig2))
sig_doub<-sig2+sig1
sig_sing<-sig[,!c(new_bc %in% union(lookup$original,lookup$replacement))]
sig_full<-cbind(sig_sing,sig_doub)

sig_full<-sig_full>0
hto2<-hto_full
shift<-c(0.3,0.2,0.2,0.2,0.25,0.1)
hto_shift<-shift*hto_full
hto2[sig_full]<-round(hto_shift[sig_full])

class(hto2)
hto2<-as.matrix(hto2)

hto_low<-SingleCellExperiment(list(counts=hto2))
altExp(sce,"HTO_old")<-altExp(sce,"HTO")
altExp(sce,"HTO")<-hto_low

sig_sce<-SingleCellExperiment(list(counts=sig_full))
altExp(sce,"sig")<-sig_sce

sce$pc_rna<-1 
sce$pc_rna[sce$truth=="Doublet"]<-pc_rna

save(sce,file="sce.rdata")
