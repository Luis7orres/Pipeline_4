# SNAKEFILE PIPELINE 4.2
import os
import snakemake.io
import glob

# Load config file
configfile: "/mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/pipeline_4.2/test/config.yaml"

# Define folders
INPUT = config["general"]["input_dir"]
OUTPUT = config["general"]["output_dir"]
RESULTS = config["general"]["results"]
LOGS = config["general"]["logs_dir"]
SCRIPTS = config["general"]["scripts_dir"]

# Define reference KMERS
KMER_DATABASE=config["params"]["ref_database"]

# Define parameters
SEQUENCING=config["params"]["seq"]
WIN_SIZE=config["params"]["win_size"]
MEAN_QUAL=config["params"]["mean_qual"]

# Define samples
SAMPLES = [os.path.basename(fastq).replace('.fastq', '') for fastq in glob.glob(os.path.join(INPUT, "*.fastq"))]

# Expected final output
rule all:
    input:
        expand(os.path.join(OUTPUT, "vsearch", "chimera_deletion", "{sample}_chimera_del.fastq"), sample=SAMPLES)

# 1- DEMULTIPLEXING WITH CUTADAPT
rule cutadapt:
    input:
        fastq=os.path.join(INPUT, "{sample}.fastq") 
    output:
        os.path.join(OUTPUT, "cutadapt", "{sample}_trimmed.fastq")
    shell:
        "bash {SCRIPTS}/1-cutadapt.sh {input.fastq} {output}"

# 2 - QUALITY CONTROL AND FILTERING
rule fastp:
    input:
        trimmed=os.path.join(OUTPUT, "cutadapt", "{sample}_trimmed.fastq")
    output:
        os.path.join(OUTPUT, "fastp", "{sample}_filt.fastq")
    params:
        SEQUENCING,
        WIN_SIZE,
        MEAN_QUAL,
        html_report=os.path.join(OUTPUT, "fastp", "reports", "{sample}_fastp_report.html"),
        json_report=os.path.join(OUTPUT, "fastp", "reports", "{sample}_fastp_report.json")
    shell:
        "bash {SCRIPTS}/2-fastp.sh {input.trimmed} {output} {params}"

# 3 - CHIMERA DELETION AND DEREPLICATION
rule vsearch:
    input:
        filtered=os.path.join(OUTPUT, "fastp", "{sample}_filt.fastq")
    output:
        nonchimeric=os.path.join(OUTPUT, "vsearch", "chimera_deletion", "{sample}_chimera_del.fastq")
    shell:
        "bash {SCRIPTS}/3-vsearch.sh {input.filtered} {output.nonchimeric}"

# # 4 - MAPPING
# rule mapping:
#     input:
#         vsearch=os.path.join(OUTPUT, "vsearch", "{sample}_vsearch.fastq")
#     output:
#         os.path.join(OUTPUT, "minimap2", "{sample}_cropped.fastq")
#     params:
#         SEQUENCING, 
#         KMER_DATABASE
#     shell:
#         "bash {SCRIPTS}/4-mapping.sh {input.vsearch} {output} {params}"

# # 5 - ASSEMBLY
# rule assembly:
#     input:
#         minimap2=os.path.join(OUTPUT, "minimap2", "{sample}_cropped.fastq")
#     output:
#         os.path.join(OUTPUT, "spades", "{sample}_spades.fasta")
#     params:
#         SEQUENCING,
#         KMER_DATABASE
#     shell:
#         "bash {SCRIPTS}/5-assembly.sh {input.minimap2} {output} {params}"
