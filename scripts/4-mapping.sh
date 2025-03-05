#!/bin/bash
# filepath: /mnt/lustre/scratch/nlsas/home/uvi/bg/sbg/pipelines/Pipeline_4/scripts/4-mapping.sh

# 1. Parse arguments
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <input_fastq> <output_sam> <output_bam> <output_bed> <output_fastq> <output_report> <seq_tech> <reference>"
    exit 1
fi

INPUT=$1
OUTPUT_SAM=$2
OUTPUT_BAM=$3
OUTPUT_BED=$4
OUTPUT_FASTQ=$5
OUTPUT_REPORT=$6
SEQ_TECH=$7
REFERENCE=$8

# 2. Load modules
module load cesga/2020  
module load gcccore/system
module load k8/0.2.5
module load minimap2/2.28
module load htslib/1.19
module load samtools/1.19
module load bedtools/2.31.0

# 3. Create directories for outputs
mkdir -p "$(dirname "$OUTPUT_SAM")" "$(dirname "$OUTPUT_BAM")" "$(dirname "$OUTPUT_BED")" "$(dirname "$OUTPUT_FASTQ")" "$(dirname "$OUTPUT_REPORT")"

# 4. Validate inputs
if [ ! -s "$INPUT" ]; then
    echo "Error: Input file empty or missing: $INPUT"
    exit 1
fi

if [ ! -s "$REFERENCE" ]; then
    echo "Error: Reference file empty or missing: $REFERENCE"
    exit 1
fi

# 5. Run minimap2 based on sequencing technology
echo "Running Minimap2 for $SEQ_TECH sequencing"
if [ "$SEQ_TECH" == "nanopore" ]; then
    minimap2 -ax map-ont -k 15 -w 5 -A2 -B4 -O4,24 -E2,1 "$REFERENCE" "$INPUT" > "$OUTPUT_SAM"
elif [ "$SEQ_TECH" == "illumina" ]; then
    minimap2 -ax sr "$REFERENCE" "$INPUT" > "$OUTPUT_SAM"
else
    echo "Error: Unsupported sequencing technology: $SEQ_TECH"
    exit 1
fi

# 6. Process SAM/BAM files
samtools view -S -b "$OUTPUT_SAM" | samtools sort -o "$OUTPUT_BAM"
samtools index "$OUTPUT_BAM"

# 7. Check for valid alignments
if [ "$(samtools view -c "$OUTPUT_BAM")" -eq 0 ]; then
    echo "Warning: No valid alignments found"
    touch "$OUTPUT_BED" "$OUTPUT_FASTQ" "$OUTPUT_REPORT"
    exit 0
fi

# 8. Generate outputs
bedtools bamtobed -i "$OUTPUT_BAM" > "$OUTPUT_BED"
samtools fastq "$OUTPUT_BAM" > "$OUTPUT_FASTQ"

# 9. Generate TSV mapping report
echo -e "ReadID\tReference\tMappingQuality" > "$OUTPUT_REPORT"
samtools view "$OUTPUT_BAM" | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $5}' >> "$OUTPUT_REPORT"
echo "Mapping report generated: $OUTPUT_REPORT"

echo "Mapping completed successfully"
exit 0