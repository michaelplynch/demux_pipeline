## create simulated sce object

library(tictoc)
library(SingleCellExperiment)
library(Matrix)
library(tidyverse)
library(demuxSNP)

args <- commandArgs(trailingOnly = TRUE)
tenx_path<-args[1]#c("/data/demuxSNP/data/processed/pbmc_data_4x_out/outs/filtered_feature_bc_matrix")
barcodes_path<-args[2]#c("/data/scratch/nextflow/data_sub/barcodes_merged_p.tsv")
lookup_path<-args[3]#c("/data/scratch/nextflow/data_sub/lookup_table_doublets_pbmc_p_40pc.tsv")
key<-args[4]

matrix<-readMM(paste(tenx_path,"/matrix.mtx.gz",sep=""))
barcodes<-read.table(paste(tenx_path,"/barcodes.tsv.gz",sep=""))
features<-read.table(paste(tenx_path, "/features.tsv.gz",sep=""))

rownames(matrix)<-features$V2
colnames(matrix)<-barcodes$V1

rna<-matrix[features$V3=="Gene",]
hto<-matrix[features$V3=="Multiplexing",]

sub_barcodes<-read.table(barcodes_path)
lookup<-read_tsv(lookup_path)

rna<-rna[,substr(colnames(rna),start=1,stop = 16) %in% substr(sub_barcodes$V1, start=1,stop=16)]
hto<-hto[,substr(colnames(hto),start=1,stop = 16) %in% substr(sub_barcodes$V1, start=1,stop=16)]

new_bc<-sub_barcodes$V1[match(substr(colnames(rna),1,16),substr(sub_barcodes$V1,1,16))]
rna1<-rna[,new_bc %in% lookup$original]
rna2<-rna[,new_bc %in% lookup$replacement]
print(dim(rna1));print(dim(rna2))
rna_doub<-rna2+rna1
rna_sing<-rna[,!c(new_bc %in% union(lookup$original,lookup$replacement))]
rna_full<-cbind(rna_sing,rna_doub)

sce<-SingleCellExperiment(list(counts=rna_full))
mainExpName(sce)<-"RNA"
head(colnames(sce))

sce$new_bc<-sub_barcodes$V1[match(substr(colnames(sce),1,16),substr(sub_barcodes$V1,1,16))]

sce$truth<-substr(sce$new_bc,start=18,stop=19)
sce$truth[sce$new_bc %in% lookup$replacement]<-"Doublet"
table(sce$truth)

source('/home/m.lynch/Github/demuxSNP-paper-figures/R/hto_sim.R')

set.seed(2)
sce
ngroup=6
z<-matrix(0,nrow=6,ncol=dim(sce)[2])
for (i in seq_len(ngroup)) {
  names<-c("K1","K2","K3","K4","K5","K6")
  n<-names[i]
  z[i,sce$truth==n]<-1

}

nsinglet=as.vector(table(sce$truth[sce$truth!="Doublet"]))[seq_len(ngroup)]
ndoub=sum(sce$truth=="Doublet")
doub<-matrix(0,nrow=ngroup,ncol=ndoub)
prob=proportions(nsinglet)
ind<-replicate(n=ndoub,sample(seq_len(ngroup),size=2,prob=prob))
for (j in seq_len(ndoub)) {
  doub[ind[,j],j]<-1
}
z[,sce$truth=="Doublet"]<-doub

#size_sig=rep(10,ngroup)
#size_bg=rep(7,ngroup)
#mu_sig=c(400,500,200,400,400,400)
#mu_bg=c(90,100,60,80,60,60)

size_sig=rep(10,6)
size_bg=rep(7,6)
mu_sig=c(500,500,400,500,200,400)
mu_bg=c(110,120,90,100,60,80)

counts<-draw_counts(size_sig = size_sig, size_bg = size_bg, mu_sig = mu_sig, mu_bg = mu_bg, mat = z)
stopifnot(!is.na(counts))
colnames(counts)<-colnames(sce)
altExp(sce,"HTO")<-SingleCellExperiment(assays=list(counts=counts))
save(sce,file="sce.rdata")
