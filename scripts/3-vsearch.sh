#!/bin/bash 

# Establish input and output routes
# INPUT ARGUMENTS
INPUT_FASTA=$1 
OUTPUT_NONCHIMERIC=$2 

# Load modules
module load cesga/2020
module load gcccore/system
module load vsearch/2.17.1

# CHIMERA DELETION 
echo "Detecting and deleting chimeras"
vsearch --uchime_denovo "$INPUT_FASTA" --nonchimeras "$OUTPUT_NONCHIMERIC"

echo "Chimera deletion complete. Exiting vsearch."
