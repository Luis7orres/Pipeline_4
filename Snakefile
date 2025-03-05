# SNAKEFILE PIPELINE 4.2
import os
import snakemake.io
import glob

# Load config file
configfile: "/mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/Pipeline_4/config.yaml"

# Define folders
INPUT = config["general"]["input_dir"]
OUTPUT = config["general"]["output_dir"]
RESULTS = config["general"]["results"]
LOGS = config["general"]["logs_dir"]
SCRIPTS = config["general"]["scripts_dir"]

# Define reference KMERS
KMER_DATABASE = config["params"]["ref_database"]

# Define parameters
SEQUENCING = config["params"]["seq"]
WIN_SIZE = config["params"]["win_size"]
MEAN_QUAL = config["params"]["mean_qual"]

# Define samples
SAMPLES = [os.path.basename(fastq).replace('.fastq', '') for fastq in glob.glob(os.path.join(INPUT, "*.fastq"))]

# Expected final output
rule all:
    input:
        expand(os.path.join(OUTPUT, "1.cutadapt", "{sample}_trimmed.fastq"), sample=SAMPLES),
        expand(os.path.join(OUTPUT, "2.fastp", "{sample}_filt.fastq"), sample=SAMPLES),
        expand(os.path.join(OUTPUT, "3.vsearch", "{sample}_chimera_del.fastq"), sample=SAMPLES),
        expand(os.path.join(OUTPUT, "4.minimap2", "{sample}_mapping_report.tsv"), sample=SAMPLES),
        os.path.join(RESULTS, "kmer_counts", "kmer_counts_matrix.tsv"),
        os.path.join(RESULTS, "bacteria_matches", "bacteria_counts_matrix.tsv")

# 1- DEMULTIPLEXING WITH CUTADAPT
rule cutadapt:
    input:
        fastq = os.path.join(INPUT, "{sample}.fastq")
    output:
        os.path.join(OUTPUT, "1.cutadapt", "{sample}_trimmed.fastq")
    shell:
        "bash {SCRIPTS}/1-cutadapt.sh {input.fastq} {output}"

# 2 - QUALITY CONTROL AND FILTERING
rule fastp:
    input:
        trimmed = os.path.join(OUTPUT, "1.cutadapt", "{sample}_trimmed.fastq")
    output:
        os.path.join(OUTPUT, "2.fastp", "{sample}_filt.fastq")
    params:
        SEQUENCING,
        WIN_SIZE,
        MEAN_QUAL,
        html_report = os.path.join(OUTPUT, "2.fastp", "reports", "{sample}_fastp_report.html"),
        json_report = os.path.join(OUTPUT, "2.fastp", "reports", "{sample}_fastp_report.json")
    shell:
        "bash {SCRIPTS}/2-fastp.sh {input.trimmed} {output} {params}"

# 3 - CHIMERA DELETION AND DEREPLICATION
rule vsearch:
    input:
        filtered = os.path.join(OUTPUT, "2.fastp", "{sample}_filt.fastq")
    output:
        nonchimeric = os.path.join(OUTPUT, "3.vsearch", "{sample}_chimera_del.fastq")
    shell:
        "bash {SCRIPTS}/3-vsearch.sh {input.filtered} {output.nonchimeric}"

rule mapping:
    input:
        vsearch = os.path.join(OUTPUT, "3.vsearch", "{sample}_chimera_del.fastq")
    output:
        sam = temp(os.path.join(OUTPUT, "4.minimap2", "sam_files", "{sample}.sam")),
        bam = os.path.join(OUTPUT, "4.minimap2", "bam_files", "{sample}.bam"),
        bed = os.path.join(OUTPUT, "4.minimap2", "bed_files", "{sample}.bed"),
        final_fastq = os.path.join(OUTPUT, "4.minimap2", "reports", "{sample}_cropped.fastq"),
        report = os.path.join(OUTPUT, "4.minimap2", "{sample}_mapping_report.tsv")
    params:
        sequencing = SEQUENCING,
        reference = KMER_DATABASE
    log:
        os.path.join(LOGS, "4.mapping", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.sam}) $(dirname {output.bam}) $(dirname {output.bed}) $(dirname {output.final_fastq})
        bash {SCRIPTS}/4-mapping.sh {input.vsearch} {output.sam} {output.bam} {output.bed} {output.final_fastq} {output.report} {params.sequencing} {params.reference} 2> {log}
        """

rule table_counts:
    input:
        reports = expand(os.path.join(OUTPUT, "4.minimap2", "{sample}_mapping_report.tsv"), sample=SAMPLES)
    output:
        counts_table = os.path.join(RESULTS, "kmer_counts", "kmer_counts_matrix.tsv"),
        bacteria_table = os.path.join(RESULTS, "bacteria_matches", "bacteria_counts_matrix.tsv")
    params:
        database = KMER_DATABASE
    script:
        os.path.join(SCRIPTS, "5-counts_table.py")
