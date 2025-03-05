# Pipeline_4: Sequence Processing & K-mer Analysis

## Overview

A Snakemake pipeline for processing sequencing data, optimized for Nanopore reads and designed for pathogen detection using the **pathogen detection kit based on probes from Science and Business SL**. This updated pipeline incorporates demultiplexing and adapter trimming, quality filtering, chimera removal, and reference mapping to accurately identify pathogens. The final step generates count matrices for k-mer and bacterial profiling.

## Table of Contents

- [Pipeline Structure](#pipeline-structure)
- [Requirements](#requirements)
- [Workflow](#workflow)
- [Rules Description](#rules-description)
- [Configuration](#configuration)
- [Config File](#config-file)
- [Usage](#usage)

## Pipeline Structure


```mermaid
flowchart TD
    A[Input FASTQ Files]
    B[1.Demultiplexing & <br>Adapter Trimming<br><i>Cutadapt</i>]:::purple
    B_out[Trimmed FASTQ:<br>sample_trimmed.fastq]
    C[2.Quality Control & Filtering<br><i>Fastp</i>]:::purple
    C_out[Filtered FASTQ:<br>sample_filt.fastq]
    D[3.Chimera Removal<br><i>Vsearch</i>]:::purple
    D_out[Chimera Removed FASTQ:<br>sample_chimera_del.fastq]
    E[4.Reference Mapping<br><i>Minimap2</i>]:::purple
    E_out[Mapping Report:<br>sample_mapping_report.tsv]:::green
    F[5.Generate K-mer Count Matrix]:::purple
    G[5.Generate Bacterial Count Matrix]:::purple
    H[K-mer Counts Matrix:<br>kmer_counts_matrix.tsv]:::green
    I[Bacteria Counts Matrix:<br>bacteria_counts_matrix.tsv]:::green
    J[Pending:<br>Epi2Me Reporting]:::white

    A --> B
    B --> B_out
    B_out --> C
    C --> C_out
    C_out --> D
    D --> D_out
    D_out --> E
    E --> E_out
    E_out --> F
    E_out --> G
    F --> H
    G --> I
    E --> J

    %% Define styles
    classDef purple fill:#734f9a,stroke:#333,stroke-width:2px,rx:5px,ry:5px,color:#fff
    classDef green fill:#3f6d4e,stroke:#333,stroke-width:2px,rx:5px,ry:5px,color:#fff
    classDef white fill:#fff,stroke:#333,stroke-width:2px,rx:5px,ry:5px,color:#000

    %% Apply styles
    class J white
    class E_out,H,I green
    class B,C,D,E,F,G purple
```

The pipeline consists of interconnected rules that handle different stages of the data processing workflow:

1. **Demultiplexing**: Using Cutadapt to trim and demultiplex the reads.
2. **Quality Control and Filtering**: Applying Fastp for quality control, filtering, and report generation.
3. **Chimera Detection and Dereplication**: Utilizing VSEARCH to detect chimeras and dereplicate sequences.
4. **Reference Mapping**: Using Minimap2 to map processed reads against a reference database, generating alignment files (SAM/BAM/BED) and a mapping report.
5. **Count Matrix Generation**: Generating two count matrices:
- *K-mer Count Matrix*: Quantifying k-mer matches from the mapping results.
- *Bacterial Count Matrix*: Tracking bacterial reference hits for taxonomic profiling.
6. **Pending: Epi2Me Reporting**: Future integration of Epi2Me for generating comprehensive analysis reports (under development).

## Requirements

- [Python 3.12.6](https://www.python.org/downloads/release/python-3126/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [Fastp](https://github.com/OpenGene/fastp)
- [VSEARCH](https://github.com/torognes/vsearch)
- [Minimap2](https://github.com/lh3/minimap2)
- [Samtools](http://www.htslib.org/)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Epi2Me](https://github.com/epi2me-labs)


## Workflow

The workflow is divided into several stages, each implemented as a rule in the Snakemake pipeline. The rules are executed in sequence, with the output of one rule serving as the input for the next.

## Rules Description

### 1. Cutadapt

- **Function**: Trim adapters from reads and demultiplex the data.
- **Input**: Raw FASTQ files.
- **Output**: Trimmed FASTQ files (fastq1, fastq2, fastq3).

---

### 2. Fastp

- **Function**: Perform quality control on the trimmed reads, including filtering based on quality scores and generating reports.
- **Input**: Trimmed FASTQ files from Cutadapt.
- **Output**: 
  - Filtered FASTQ files (filtered_fastq).
  - Quality reports (HTML and JSON).

---

### 3. VSEARCH

- **Function**: Detect chimeras in the filtered sequences and dereplicate them.
- **Input**: Filtered FASTQ files from Fastp.
- **Output**: 
  - Non-chimeric FASTQ files (fastq_no_chimeras).
  - Dereplicated sequences.

---

### 4. Minimap2

- **Function**: Align the non-chimeric reads and prepare for consensus generation.
- **Input**: The fastq file with no chimeras.
- **Output**: 
  - Aligned FASTQ files (fastq_aligned_trimmed).
  - Read position information (read_positions.bed).

## Configuration

The pipeline's parameters can be configured in the `config.yaml` file. Key parameters include:

### Config File

The configuration file requires the following structure to define input and output directories, logging, and result paths, as well as analysis parameters such as sequencing type and k-mer reference database.

```bash
general:
  input_dir: /your/data/directoy
  output_dir: /desired/output/directory
  scripts_dir: /pipeline/scripts/directory
  logs_dir: /desired/logs/directory
  results: /desired/results/directory

params:
  seq: "single-end"  # or "paired-end"
  win_size: 4
  mean_qual: 7
  ref_database: /pathogen/kmer/database
```
### Key Parameters

- **input_dir**: Directory containing the input sequencing data.
- **output_dir**: Directory where the output files will be generated.
- **scripts_dir**: Directory containing the scripts used in the pipeline.
- **logs_dir**: Directory for storing log files generated during the pipeline execution.
- **results**: Directory for storing the final results of the analysis.
- **seq**: Type of sequencing data being processed (e.g., `single-end` or `paired-end`).
- **win_size**: Window size used for quality filtering.
- **mean_qual**: Minimum mean quality score for reads to be retained.
- **ref_database**: Path to the reference database for k-mers.

## Usage

To run the pipeline, navigate to the directory containing the `Snakefile` and execute:

```bash
snakemake --cores <number_of_cores>
```
Replace <number_of_cores> with the desired number of processing threads.

## THIS PIPELINE IS OPTIMIZED FOR FINISTERRAE III, WE PLAN TO CONTAINEREIZE IT WITH DOCKER TO ENSURE REPRODUCIBILITY ACROSS DIFFERENT ENVIROMENTS
