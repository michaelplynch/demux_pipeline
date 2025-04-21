args <- commandArgs(trailingOnly = TRUE)
barcodes<-Sys.glob(args[1])
print(barcodes)

library(tools)
library(readr)
library(BiocGenerics)
for (i in seq_along(barcodes)) {
  barcode<-barcodes[i]
  if (grepl('.csv',barcode)) {
    bcs<-read_csv(barcode,col_names=c('genome','barcode'))
    bcs<-bcs[,'barcode']
  } else if (grepl('.tsv',barcode)) {
    bcs<-read_table(barcode, col_names = "barcode")
  } else {print('unknown file extension')}
file_path<-c(paste0(file_path_sans_ext(barcode),'.tsv'))
  print(file_path)
  write_tsv(bcs,file=file_path,col_names = F)
}
