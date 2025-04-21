## script to carry out demuxSNP preprocessing steps

## inputs
args <- commandArgs(trailingOnly = TRUE)
vcf_path<-args[1]#c("/data/scratch/nextflow/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf")
sce_path<-args[2]#c("/home/m.lynch/Github/demuxSNP-paper-figures/sce_p.rds")
key<-args[3]
doublets<-args[4]
seed<-args[5]
n_genes<-args[6]
##

library(demuxSNP)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
# load data
load(sce_path)

# find top genes
top_genes<-common_genes(sce,n=n_genes)

# subset to 5'
# credit to Peter Hickey https://github.com/lmweber/snp-dmx-cancer/blob/master/filter_vcf/filter_vcf_mod.R
# library(AnnotationHub)
# library(ensembldb)
# ah <- AnnotationHub()
# EnsDb.Hsapiens.v94 <- ah[["AH64923"]]
# five_utrs <- fiveUTRsByTranscript(EnsDb.Hsapiens.v94)
# five_utrs <- keepSeqlevels(five_utrs, c(1:22, "X"), pruning.mode = "coarse")
# reduced_five_utrs <- reduce(unlist(five_utrs))

# subset to genes expressed
# vcf<-readVcf(vcf_path,genome="GRCh38",param = ScanVcfParam(
#   fixed = names(fixed(scanVcfHeader(vcf_path))),
#   info = names(info(scanVcfHeader(vcf_path))),
#   geno = names(geno(scanVcfHeader(vcf_path))),
#   which = reduced_five_utrs))
vcf<-readVcf(vcf_path,genome="GRCh38")
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
seqlevelsStyle(ensdb)<-"UCSC"
vcf_sub<-subset_vcf(vcf, top_genes, ensdb)

# print new vcf
dir<-paste("demuxSNP_",key,"_",doublets,"_",seed,sep="")
dir.create(dir)
writeVcf(vcf_sub,filename=paste(dir,"/vcf_sub.vcf",sep=""))
getwd()
