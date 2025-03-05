#!/usr/bin/env python3
import subprocess
import sys
print(sys.version)

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

try:
    import pandas as pd
except ImportError:
    install("pandas")
    import pandas as pd

from pathlib import Path

def extract_species(header):
    """
    Extract the bacterial species name from the FASTA header.
    This function:
      - Splits the header (without the '>') by underscores, removing empty tokens.
      - Always takes the first two tokens (assuming genus and species).
      - Then iterates over subsequent tokens, adding them if they are considered part of the species name.
      - Stops adding tokens if a metadata marker is encountered:
          * The token is "GCF" or "kmer".
          * If the third token is numeric (e.g., "770") or a strain marker ("str" or "st").
      - This approach allows keeping cases like "escherichia coli o25" or "clostridium botulinum cepa 2".
    """
    tokens = [t for t in header.split("_") if t]  # Remove empty tokens
    species_tokens = []
    for i, token in enumerate(tokens):
        # Stop if a metadata marker is encountered
        if token in ("GCF", "kmer"):
            break
        # Always include the first two tokens (genus and species)
        if i < 2:
            species_tokens.append(token)
        else:
            # For the third token: if it's purely numeric or a strain marker, stop adding further tokens.
            if i == 2 and (token.isdigit() or token.lower() in ("str", "st")):
                break
            species_tokens.append(token)
    return " ".join(species_tokens)

def main(input_reports, reference_db, output_file, bacteria_output):
    # Load the kmer reference FASTA and map each sequence to its header (without ">")
    kmer_mapping = {}
    with open(reference_db) as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]  # Remove ">"
            else:
                if header is not None:
                    # Assume the sequence is contained in a single line.
                    kmer_mapping[line] = header
                    header = None

    # Count dictionary: keys are full kmer headers (names) and values are dictionaries with sample counts.
    counts_dict = {name: {} for name in kmer_mapping.values()}

    # Process every report file
    for report in input_reports:
        # Get sample name from the report file name (assumes sample is before the first underscore)
        sample = Path(report).name.split('_')[0]

        try:
            # Read the report ignoring the header; assume column 1 contains the sequence.
            df = pd.read_csv(report, sep='\t', usecols=[1], skiprows=1, header=None)
            counts = df[1].value_counts()

            # Update counts: for each kmer (sequence) use its corresponding header and assign its count
            for seq, name in kmer_mapping.items():
                counts_dict[name][sample] = counts.get(seq, 0)
        except pd.errors.EmptyDataError:
            # If the file is empty, assign a count of 0 for each kmer.
            for name in counts_dict:
                counts_dict[name][sample] = 0
    
    # Convert the dictionary to a DataFrame and fill missing values with 0.
    result_df = pd.DataFrame.from_dict(counts_dict, orient='index')
    result_df.fillna(0, inplace=True)
    result_df = result_df.astype(int)

    # Sort by kmer names and sample names.
    result_df.sort_index(inplace=True)
    result_df = result_df[sorted(result_df.columns)]

    # Save the kmer count table.
    result_df.to_csv(output_file, sep='\t', header=True, index=True)

    # Generate an aggregated count table by bacteria:
    # Extract the species from each header and group by species, summing the counts.
    df_species = result_df.copy()
    df_species.index = [extract_species(idx) for idx in df_species.index]
    bacteria_aggregated = df_species.groupby(df_species.index).sum()
    bacteria_aggregated.sort_index(inplace=True)

    # Save the aggregated bacteria count table.
    bacteria_aggregated.to_csv(bacteria_output, sep='\t', header=True, index=True)

if __name__ == "__main__":
    main(
        input_reports=snakemake.input.reports,
        reference_db=snakemake.params.database,
        output_file=snakemake.output.counts_table,
        bacteria_output=snakemake.output.bacteria_table
    )
