# Portable genotype-free demultiplexing benchmarkign pipeline.

A portable pipeline for benchmarking genotype-free single-cell demultiplexing methods on simulated data.

The pipeline is designed to be generelisable to different datasets with arbitrary numbers of simulated mulitplexed samples.
All software as part of pipeline is run through apptainer containers to ensure reproducibility and ease of use.
The pipeline is designed to be run on a cluster with a slurm scheduler, but can be run locally with minor modifications.

The pipeline is adapted from that used to benchmark demuxSNP (DOI:) to be generalizable to other methods, datasets and research questions.
In the absence of hashing simulations techniques, the pipeline only covers genotype-based methods, not hybrid.

The pipeline can be tested using an [example dataset](https://www.10xgenomics.com/datasets/10-k-1-1-mixture-of-raji-and-jurkat-cells-multiplexed-2-cm-os-3-1-standard-6-0-0) from 10X Genomcis consisting of cells from two cells lines.

## Overall workflow

1. Simulate doublets  
Scripts/templates leveraged from Weber et al. (DOI: https://doi.org/10.1093/gigascience/giab062)
- Start with an abritrary number of demultiplexed/individual bams
- Barcode suffix is replace with sample key e.g. K1, K2 etc.
- bam files are merged.
- Lookup file generated with barcodes randomly selected to satisfy required proportion of doublets.
- Barcodes in bam renamed according to lookup to simulate doublets.  
2. Benchmark methods  
Tests demuxSNP (hybrid), HTOreader (hybrid) and souporcell.
- Experiments 1: Vary doublet rate
- Experiment 2: Vary SNP subsetting

## Inputs

Most inputs are specified in `nextflow.config`.
- `container__souporcell`: Path to the souporcell Apptainer image, ideally at the top level of the project.
- `bam_path`: Path to an arbitrary number of (demultiplexed) BAM files.
- `barcodes_path`: Path to the corresponding barcodes files (.csv).
- `common_variants`: Common variants, e.g., from the 1K Genome Project.
- `ref`: Path to the reference genome, ideally in the `data/input` directory.

Doublet simulation parameters are specified in params.csv
The workflow caters for subsampling (also specified in params.csv).

## Outputs

Output folder for each method applied to eachsimulated scenario (e.g. seed, key)

## Known issues

- Apptainer must be bound to the project directory (set in Nextflow.config).
- The pipeline runs by Slurm as default. Each process is submitted as a job once all inputs are available. On busy clusters, this may result in significant time spent in the queue, and running the pipeline within a single job may be more efficient.
- If the cluster requires apptainer to be load through modules environment, it must be loaded prior to executing the pipeline, otherwise the apptainer command is not available to the script to pull the images.