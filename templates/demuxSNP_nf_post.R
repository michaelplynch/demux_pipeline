## post vartrix demuxSNP
sink("r.log", type=c("output", "message"))
library(demuxSNP)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(ComplexHeatmap)
library(tidyverse)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
sce_path<-args[1] #c("/home/m.lynch/Github/demuxSNP-paper-figures/sce_p.rds")
snps_path<-args[2] #c("/home/m.lynch/out_matrix.mtx")
barcodes_path<-args[3] #c("/data/scratch/nextflow/data_sub/barcodes_merged_pbmc_p_40pc.tsv")
key<-args[4]
doublets<-args[5]
seed<-args[6]
souporcell<-args[7]

# load data
load(sce_path)
mat<-readMM(snps_path)
barcodes<-read_tsv(barcodes_path,col_names=FALSE)

mat_order<-as.matrix(mat[,match(substr(sce$new_bc,1,16),substr(barcodes$X1,1,16))])

sce<-add_snps(sce,mat=mat_order,thresh=0.4)
set.seed(seed)
sce<-high_conf_calls(sce, pacpt = 0.7)
sce$labels<-as.character(sce$labels)
sce$labels[sce$labels=="multiplet"]<-"Doublet"
sce$labels<-as.factor(sce$labels)
altExp(sce,"SNP")
table(sce$train,sce$truth)
table(sce$train,sce$labels)
predict = sce$labels=="uncertain" | sce$labels=="negative"
set.seed(seed)
sce<-reassign_balanced(sce,k=25,d_prop=0.5,predict_cells = predict)
sce<-reassign_jaccard(sce,k=25,d=200,predict_cells = predict)
sce<-reassign(sce,k=25,d=200,predict_cells = predict)

sce$predict<-TRUE
sce$knn_centroid_all<-reassign_centroid(sce,predict_cells=sce$predict)
sce$predict<-!sce$train
sce$knn_centroid_nottrain<-reassign_centroid(sce,predict_cells=sce$predict)
sce$predict<-sce$labels=="uncertain"|sce$labels=="negative"
sce$knn_centroid_neg_uncertain<-reassign_centroid(sce,predict_cells=sce$predict)

Heatmap(counts(altExp(sce,"SNP")),
        column_split = sce$truth,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE)

library(mclust)
adjustedRandIndex(sce$knn,sce$truth)

library(HTOreader)
# load souporcell results
souporcell.res <- read.csv(paste(souporcell,'/clusters.tsv', sep = ''),sep='\t',header=TRUE,stringsAsFactors=FALSE)
souporcell.res$identity <- souporcell.res$status
souporcell.res$identity[which(souporcell.res$status == 'singlet')] <- paste0(souporcell.res$status[which(souporcell.res$status == 'singlet')], souporcell.res$assignment[which(souporcell.res$status == 'singlet')])
rownames(souporcell.res) <- souporcell.res$barcode
souporcell.res<-souporcell.res[colnames(sce),]
sce$souporcell <- souporcell.res$identity
# run HTOreader
sce2<-sce
altExp(sce2,"sig")<-NULL;altExp(sce2,"HTO_old")<-NULL;altExp(sce2,"SNP")<-NULL
seurat<-as.Seurat(sce2[unique(rownames(sce2)),],data=NULL)
set.seed(seed)
i=0;s2<-c()
while (class(s2)!="Seurat" & i<20) {
        try({
                s2 <- HTOClassification(seurat, assay = 'HTO', method = 'log')
                PlotHTO(s2, assay = 'HTO', method = 'log')
                s2 <- HybridDemultiplexing(s2, cellhashing_label = 'HTOid', genotype_label = 'souporcell', hto_names = paste0("Hashtag",1:6),c_threshold = 0)
        }
        )
  i=i+1
  print(i)
}
seurat<-s2
sce$htoreader<-seurat$HTOid
sce$htoreadermulti<-seurat$hybridID

table(sce$htoreadermulti,sce$htoreader)

result<-data.frame(barcode=colnames(sce),
        demuxSNP=sce$knn,
        demuxSNP_jacc=sce$knn_jacc,
        demuxSNP_balanced=sce$knn_balanced,
        htoreader=sce$htoreader,
        htoreadermulti=sce$htoreadermulti,
        souporcell=sce$souporcell,
        demuxSNP_centroid_all=sce$knn_centroid_all,
        demuxSNP_centroid_neguncert=sce$knn_centroid_neg_uncertain,
        demuxSNP_centroid_nottrain=sce$knn_centroid_nottrain)
head(result)

dir<-paste("demuxSNP_",key,"_",doublets,"_",seed,sep="")
dir.create(dir,recursive=TRUE)
write.table(result,file=paste(dir,"/",key,"_demuxSNP.tsv",sep=""),col.names = TRUE,quote=FALSE,row.names = FALSE)
save(sce,file=paste(dir,"/",key,"_final_sce.rdata",sep=""))
