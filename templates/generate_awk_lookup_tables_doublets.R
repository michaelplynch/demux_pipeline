#############################################################################
# Doublets simulations: generate awk lookup tables and updated barcodes files
#############################################################################

# Simulate doublets by combining some percentage of cell barcodes in the merged
# BAM file.

# This script generates:
# (i) a lookup table for the simulation scenario, containing the set of cell
# barcodes to replace and combine. This lookup table is then used in the awk
# command in the shell script for the simulation scenario.
# (ii) updated merged cell barcodes file for each scenario


# module load conda_R/4.0
# qsub -V -cwd -l mem_free=2G,h_vmem=3G,h_fsize=100G Rscript generate_awk_lookup_tables_doublets.R


# ----------------------
# Command-line arguments
# ----------------------

args <- commandArgs(trailingOnly = TRUE)
dir_runtimes <- args[1]
dir_timestamps <- args[2]
dir_out <- args[3]
barcodes <- args[4]
pc_doublets_sims <- as.numeric(args[5])
key <- args[6]
dataset_name_sims <- 'ccrcc' # args[6]


# ------------------------------------------
# Function to generate lookup tables for awk
# ------------------------------------------

# function to save lookup table for each simulation scenario

# arguments:
# prop_doublets: proportion of doublets to simulate
# dataset_name: dataset name for output files
# file_barcodes_merged: tsv file containing merged list of cell barcodes from Cell Ranger
f_sim_doublets <- function(pc_doublets, dataset_name, file_barcodes_merged,key) {
  
  library(tidyverse)
  library(magrittr)
  
  # corrected proportion of doublets to simulate, i.e. number of cells to
  # replace after taking into account reduction in final number of cells
  pc_doublets_corrected <- (pc_doublets) / (100 + pc_doublets)
  print(paste('corrected # doublets=', pc_doublets_corrected))
  
  
  # ---------------------
  # Generate lookup table
  # ---------------------
  
  # load merged list of cell barcodes
  df_barcodes <- read_table(file_barcodes_merged, col_names = "barcode")
  
  # number of cells
  n_cells <- nrow(df_barcodes)
  # number of cells to combine with other cells as doublets
  n_doublets <- round(pc_doublets_corrected * n_cells)
  
  print(paste('# cells=', n_cells))
  print(paste('# doublets=', n_doublets))
  
  # select random sets of (i) cells to replace and (ii) replacements
  set.seed(123)
  ix_original <- sample(seq_len(n_cells), n_doublets)
  ix_replacement <- sample(setdiff(seq_len(n_cells), ix_original), n_doublets)
  
  # generate lookup table / replacement table for cell barcodes
  df_lookup <- tibble(
    original = df_barcodes$barcode[ix_original], 
    replacement = df_barcodes$barcode[ix_replacement]
  )
  
  # create output directory
  #dir_tmp <- paste0(dir_out, "/", dataset_name, "/doublets_sims/", prop_doublets * 100, "pc")
  #system(paste0("mkdir -p ", dir_tmp))
  outpath = dir_out
  # save lookup table
  write_tsv(df_lookup, paste0("lookup_table_doublets_", dataset_name, "_", key, "_", pc_doublets, "pc.tsv"))
  
  
  # ------------------------------
  # Generate updated barcodes file
  # ------------------------------
  
  # generate updated list of cell barcodes by removing replaced cells
  
  barcodes_merged_new <- df_barcodes$barcode[!(df_barcodes$barcode %in% df_lookup$original)]
  print(length(barcodes_merged_new))
  stopifnot(all(barcodes_merged_new %in% df_barcodes$barcode))
  stopifnot(length(barcodes_merged_new) == n_cells - n_doublets)
  
  # save barcodes file
  write_tsv(tibble(barcodes_merged_new), paste0("barcodes_merged_", dataset_name, "_", key, "_", pc_doublets, "pc.tsv"), col_names = FALSE)
  
}


# ----------------------------------------------------------------
# Save lookup table and barcodes file for each simulation scenario
# ----------------------------------------------------------------

# run function to save lookup table and updated barcodes file for each
# simulation scenario

for (pc_doublets in pc_doublets_sims) {
  for (dataset_name in dataset_name_sims) {
    runtime <- system.time({
      dir_out<-dir_out
      file_barcodes_merged <- file.path(barcodes)
      f_sim_doublets(pc_doublets, dataset_name, file_barcodes_merged, key)
    })
  }
}

