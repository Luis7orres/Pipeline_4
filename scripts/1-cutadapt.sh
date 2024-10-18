#!/bin/bash 

# Establish input and output routes
INPUT_FASTQ=$1
OUTPUT_TRIMMED=$2

# Establish the route to the fasta file containign the adapters
# ADAPTERS_FASTA=$3

#Ensure output directory exists
mkdir -p "$(dirname "$OUTPUT_TRIMMED")"

# Clean modules
module purge

# Load modules 
module load cesga/2020
module load gcccore/system
module load cutadapt/3.5

# DEMULTIPLEXING WITH CUTADAPT
cutadapt \
    -o "$OUTPUT_TRIMMED" \
    "$INPUT_FASTQ" \
    # -g file:"$ADAPTERS_FASTA" \

if [ $? -ne 0 ]; then 
    echo "Cutadapt failed, exiting pipeline 4"
    exit 1 
fi

echo "Demultiplexing finished, results located in $OUTPUT_TRIMMED"