#!/bin/bash

# Establish input and ouptut routes
INPUT_TRIMMED=$1 
OUTPUT_FASTP=$2
HTML_REPORT=$6
JSON_REPORT=$7

# Establish parameters
SEQUENCING=$3
WINDOW_SIZE=$4  
MEAN_QUAL=$5


# LOAD MODULES
module load cesga/2020  
module load gcccore/system
module load fastp/0.22.0

# FASTP FOR SINGLE END SEQUENCING

if [ "$SEQUENCING" == "single-end" ]; then
    echo "Running fastp for single end sequencing"

    # Create report folder (if necessary)
    REPORTS_DIR=$(dirname "$HTML_REPORT")
    mkdir -p "${REPORTS_DIR}" 

    # Run fastp with specified paths for the reports
    fastp -i "$INPUT_TRIMMED" -o "$OUTPUT_FASTP" \
      -W "$WINDOW_SIZE" -M "$MEAN_QUAL" \
      --cut_front --cut_tail \
      --html "$HTML_REPORT" \
      --json "$JSON_REPORT"

    echo "FASTP analysis for single-end sequencing complete"

elif [ "$SEQUENCING" == "paired-end" ]; then
    echo "Running fastp for paired-end sequencing"

    # Automatically derive R2 from R1 filename (replace _R1 with _R2)
    INPUT_TRIMMED_R2="${INPUT_TRIMMED/_R1/_R2}"

    # Output for R2 (replace _R1 with _R2 in output)
    OUTPUT_FASTP_R2="${OUTPUT_FASTP/_R1/_R2}"

    # Run fastp
    fastp -i "$INPUT_TRIMMED" -I "$INPUT_TRIMMED_R2" \
          -o "$OUTPUT_FASTP" -O "$OUTPUT_FASTP_R2" \
          -W $WINDOW_SIZE -M $MEAN_QUAL -l $MIN_LEN \
          --cut_front --cut_tail 

    echo "FASTP analysis for paired-end sequencing complete"

else
    echo "Sequencing method not recognized: $SEQUENCING"
fi

