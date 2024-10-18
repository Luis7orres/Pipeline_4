# README: Pathogen Detection Pipeline

## Overview

This pipeline is designed for processing sequencing data from Illumina and Nanopore platforms, specifically for pathogen detection using the **pathogen detection kit based on probes from Science and Business**. It provides a streamlined workflow for data demultiplexing, quality control, chimera detection, and dereplication of sequencing reads.

## Table of Contents

- [Pipeline Structure](#pipeline-structure)
- [Requirements](#requirements)
- [Workflow](#workflow)
- [Rules Description](#rules-description)
- [Configuration](#configuration)
- [Config File](#config-file)
- [Usage](#usage)
- [Output](#output)
- [License](#license)

## Pipeline Structure


```mermaid
flowchart TD
    %% Define nodes
    A[Nanopore data] 
    B[Illumina paired-end data] 

    demultiplex[Cutadapt<br><i>Demultiplexing</i>]:::color1
    fastq1[fastq 1]
    fastq2[fastq 2]
    fastq3[fastq 3]

    fastq4[fastq 4]
    fastq5[fastq 5]
    fastq6[fastq 6]

    QC[FastP<br><i>Quality Control</i>]:::color2
    filtered_fastq[filtered_fastq]
    
    QC_HTML[HTML Report]
    QC_JSON[JSON Report]

    Chimera[VSearch<br><i>Chimera Deletion</i>]:::color3
    fastq_no_chimeras[fastq_no_chimeras]

    Consensus[Minimap2<br><i>Consensus Preparation</i>]:::color4
    fastq_aligned_trimmed[fastq_aligned_trimmed]
    read_positions_bed[read_positions.bed]

    %% Define styles for colored nodes
    classDef color1 fill:#ffcc00,stroke:#333,stroke-width:2px;
    classDef color2 fill:#66ccff,stroke:#333,stroke-width:2px;
    classDef color3 fill:#ff6699,stroke:#333,stroke-width:2px;
    classDef color4 fill:#99cc66,stroke:#333,stroke-width:2px;

    %% Define links
    A --> demultiplex
    demultiplex --> fastq1
    demultiplex --> fastq2
    demultiplex --> fastq3

    B --> fastq4
    B --> fastq5
    B --> fastq6

    fastq1 --> QC
    fastq2 --> QC
    fastq3 --> QC
    fastq4 --> QC
    fastq5 --> QC
    fastq6 --> QC

    QC --> filtered_fastq

    filtered_fastq --> Chimera
    Chimera --> fastq_no_chimeras

    fastq_no_chimeras --> Consensus
    Consensus --> fastq_aligned_trimmed
    Consensus --> read_positions_bed

    %% Position the reports
    QC --> QC_HTML
    QC --> QC_JSON
```

The pipeline consists of several interconnected rules that handle different stages of the data processing workflow:

1. **Demultiplexing**: Using Cutadapt to trim and demultiplex the reads.
2. **Quality Control and Filtering**: Applying Fastp for quality control, filtering, and report generation.
3. **Chimera Detection and Dereplication**: Utilizing VSEARCH to detect chimeras and dereplicate sequences.

## Requirements

- Snakemake
- Cutadapt
- Fastp
- VSEARCH
- Access to Illumina and Nanopore sequencing data

## Workflow

The workflow is divided into several stages, each implemented as a rule in the Snakemake pipeline. The rules are executed in sequence, with the output of one rule serving as the input for the next.

## Rules Description

### 1. Cutadapt

- **Function**: Trim adapters from reads and demultiplex the data.
- **Input**: Raw FASTQ files.
- **Output**: Trimmed FASTQ files.

### 2. Fastp

- **Function**: Perform quality control on the trimmed reads, including filtering based on quality scores and generating reports.
- **Input**: Trimmed FASTQ files from Cutadapt.
- **Output**: Filtered FASTQ files and quality reports (HTML and JSON).

### 3. VSEARCH

- **Function**: Detect chimeras in the filtered sequences and dereplicate them.
- **Input**: Filtered FASTQ files from Fastp.
- **Output**: Non-chimeric FASTQ files and dereplicated sequences.

## Configuration

The pipeline's parameters can be configured in the `config.yaml` file. Key parameters include:

### Config File

The configuration file requires the following structure to define input and output directories, logging, and result paths, as well as analysis parameters such as sequencing type and k-mer reference database.

```bash
general:
  input_dir: /mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/test_pipeline4/data
  output_dir: /mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/test_pipeline4/output
  scripts_dir: /mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/Pipeline_4/scripts
  logs_dir: /mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/test_pipeline4/logs
  results: /mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/test_pipeline4/results

params:
  seq: "single-end"  # or "paired-end"
  win_size: 4
  mean_qual: 7
  ref_database: /mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/Pipeline_4/database.db
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

## Output

The processed files will be organized in the specified output directory. The structure includes directories for each processing step, with outputs named according to the sample identifiers.

- **Cutadapt Output**: Trimmed FASTQ files located in the `cutadapt` directory.
- **Fastp Output**: Filtered FASTQ files and associated reports (HTML and JSON) located in the `fastp` directory.
- **VSEARCH Output**: Non-chimeric sequences and dereplicated sequences located in the `vsearch` directory.
- ***TO BE CONTINUED***


