# RNA-Seq_Analysis

### About

This generalized pipeline will take user defined raw paired-end fastq files from illumina and run the following programs.

1) Run fastqc quality control check of each sample

2) Index the user defined transcript assembly using salmon

3) Generate pseudocounts from the .fastq files using the salmon index


## Installation

This pipeline requires fastqc and salmon.
We recommend the use of conda environments to install these packages,
following these [instructions](https://docs.conda.io/en/latest/miniconda.html) for installation of miniconda.

To create the enviornment, run the following command:

```
conda create -n rnaseq -c conda-forge -c bioconda fastqc salmon

```

To activate this enviornment:

```
conda activate rnaseq
```

## Quickstart

In order to use this pipeline, the user needs to supply a transcriptome file.
Furthermore, illumina paired-end sequencing reads should have the standard \_R1_ & \_R2_ file naming convention.

To run the pipline, run the following command

```
bash rna_seq_pipeline.sh 
-s <path/to/sequence/directory>
-o <path/to/output/directory>
-a <path/to/assembly/file>
-t <number of threads>
```
