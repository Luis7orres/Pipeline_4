#!/bin/bash
INPUT=$1
OUTPUT_MAPPING=$2
SEQ_TECH=$3

if [ $SEQ_TECH == "nanopore" ];  then
    minimap2 -ax sr "$INPUT" > aligned_reads.sam

#Mapear las lecturas contra la referencia (Illumina o Nanopore)
minimap2 -ax sr reference.fasta reads.fastq > aligned_reads.sam
minimap2 -ax sr reference.fasta reads_R1.fastq reads_R2.fastq > aligned_reads.sam

#Convertir SAM a BAM y ordenar
samtools view -bS aligned_reads.sam | samtools sort -o aligned_reads_sorted.bam

#Extraer las lecturas alineadas y las posiciones de mapeo
samtools view aligned_reads_sorted.bam | awk '{print $1, $4, $9}' > read_positions.txt

# Convertir read_positions.txt a archivo BED
awk '{print $2 "\t" $3-1 "\t" $3 }' read_positions.txt > regions.bed

#Generar un fastq con las lecturas mapeadas recortadas
bedtools intersect -a aligned_reads_sorted.bam -b regions.bed | bedtools bamtofastq -i stdin -fq recortadas.fastq