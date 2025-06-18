# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 13:51:41 2025

@author: paulo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt # Although not directly used in the functions, keeping for completeness if intended later
import argparse
import os, sys

def isBiallelic(minor_exp):
    """
    Categorizes minor allele expression into 'biallelic', 'imbalance', or 'monoallelic'.
    """
    if np.isnan(minor_exp):
        return np.nan
    elif minor_exp >= 0.40:
        return 'biallelic'
    elif minor_exp > 0.10:
        return 'imbalance'
    else:
        return 'monoallelic'

def gene_level_aggregation(ase_res, threshold = 10):
    """
    Aggregates counts from all SNPs for a given gene and computes MAF
    taking the total count into account.

    Args:
        ase_res (pd.DataFrame): Input DataFrame containing SNP-level ASE results.
        threshold (int): Minimum total counts for a SNP to be included in aggregation.

    Returns:
        pd.DataFrame: A DataFrame aggregated at the gene level with allelic expression labels.
    """
    # Apply threshold for min total counts
    # Using .copy() to avoid SettingWithCopyWarning in pandas
    ase_res_filtered = ase_res[ase_res['totalCount'] > threshold].copy()

    if ase_res_filtered.empty:
        print(f"Warning: After applying 'totalCount' threshold > {threshold}, no data remains. Returning an empty DataFrame.", file=sys.stderr)
        # Return a DataFrame with the expected columns, even if empty, to prevent errors downstream
        return pd.DataFrame(columns=[
            'gene_id', 'gene_name', 'gene_chr', 'gene_start', 'gene_end', 'gene_strand',
            'refCount', 'altCount', 'totalCount', 'minor_fq', 'variant', 'minor_fq_tot', 'allelic_exp'
        ])

    # Create variant column by concatenating relevant SNP information
    variant_label = (ase_res_filtered['contig'] + '_' +
                     ase_res_filtered['position'].astype(str) + '_' +
                     ase_res_filtered['refAllele'] + '_' +
                     ase_res_filtered['altAllele'])
    ase_res_filtered.insert(ase_res_filtered.shape[1], 'variant', variant_label)

    # Compute minor allele frequency for each SNP
    minor_fq = ase_res_filtered[['refCount','altCount']].min(axis = 1) / (ase_res_filtered['totalCount'])
    ase_res_filtered.insert(ase_res_filtered.shape[1], 'minor_fq', minor_fq)

    # Aggregate ref and alt count per gene;
    # compute mean MAF and aggregate all variants into a single column
    print("Aggregating SNP data to gene level...")
    gene_level = (
        ase_res_filtered.groupby(['gene_id', 'gene_name',
                                  'gene_chr','gene_start','gene_end', 'gene_strand'])
        .agg({
            'refCount': 'sum',
            'altCount': 'sum',
            'totalCount':'sum',
            'minor_fq': 'mean', # Mean of individual SNP minor frequencies
            'variant': lambda x: ','.join(x.astype(str)) # Concatenate all variants
        }).reset_index()
    )

    # Compute gene-level minor allele frequency based on summed counts
    # This calculation should use the aggregated refCount and altCount at the gene level
    minor_fq_tot = gene_level[['refCount','altCount']].min(axis = 1) / (gene_level['totalCount'])
    # Insert at a logical position, e.g., before allelic_exp
    gene_level.insert(gene_level.shape[1], 'minor_fq_tot', minor_fq_tot)


    # Add label according to the minor allele expression (monoallelic vs biallelic)
    # This uses the 'minor_fq' (mean SNP minor frequency) for classification
    gene_level.insert(gene_level.shape[1], 'allelic_exp', gene_level['minor_fq'].apply(isBiallelic))
    
    print("Gene level aggregation complete.")
    return gene_level

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        A command-line tool to perform gene-level aggregation of Allele-Specific Expression (ASE) results.

        This script takes an input CSV or TSV file containing SNP-level ASE data,
        filters by total read count, computes minor allele frequencies,
        aggregates data to gene level, and categorizes allelic expression.
        """,
        formatter_class=argparse.RawTextHelpFormatter # For better help formatting
    )

    # Define the input file argument
    parser.add_argument(
        'input_file',
        type=str,
        help="Path to the input file (CSV or TSV) containing SNP-level ASE results."
             "\nExpected columns:\n"
             "  - 'gene_id', 'gene_name', 'gene_chr', 'gene_start', 'gene_end', 'gene_strand'\n"
             "  - 'refCount', 'altCount', 'totalCount'\n"
             "  - 'contig', 'position', 'refAllele', 'altAllele'"
    )

    # Optional: Add an output file argument
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None, # Default is None, allowing custom logic if not provided
        help="Optional: Path to the output CSV/TSV file where the aggregated gene-level results "
             "will be saved. If not provided, the output file will be named "
             "'<input_filename>.geneLevelMAF.<original_extension>'. "
             "If output file extension is not recognized, results will be printed to console."
    )

    # Optional: Add a threshold argument
    parser.add_argument(
        '-t', '--threshold',
        type=int,
        default=10, # Default value for threshold
        help="Optional: Minimum total read count threshold for SNPs to be included "
             "in aggregation. SNPs with 'totalCount' <= threshold will be excluded. "
             "Default: 10."
    )

    args = parser.parse_args()

    # Access the arguments
    input_filepath = args.input_file
    output_filepath = args.output # This will be None if -o is not used
    threshold = args.threshold

    print(f"--- Starting Gene Aggregation ---")
    print(f"Input file: '{input_filepath}'")
    print(f"Threshold for totalCount: {threshold}")

    # 1. Validate input file and determine separator
    if not os.path.exists(input_filepath):
        print(f"Error: Input file not found at '{input_filepath}'", file=sys.stderr)
        sys.exit(1)

    # Get file components to handle default output filename and separator
    input_dir = os.path.dirname(input_filepath)
    input_basename = os.path.basename(input_filepath)
    filename_without_ext, file_extension = os.path.splitext(input_basename)
    file_extension = file_extension.lower() # Ensure lowercase for consistent checks

    if file_extension == '.csv':
        separator = ','
    elif file_extension in ('.tab', '.tsv', '.txt'): # .txt files are often tab-separated
        separator = '\t'
    else:
        print(f"Warning: Input file format '{file_extension}' not explicitly supported (.csv, .tsv, .txt). "
              "Attempting to read as CSV by default. Output will be printed to console unless explicit output file path provided.", file=sys.stderr)
        separator = ',' # Default to CSV separator if unsure

    # 2. Determine output file path (defaulting if not provided)
    output_to_console = False
    if output_filepath is None: # If user did not provide an output file path
        if file_extension in ('.tab', '.csv', '.tsv', '.txt'):
            default_output_filename = f"{filename_without_ext}.geneLevelMAF{file_extension}"
            output_filepath = os.path.join(input_dir, default_output_filename)
            print(f"No output file specified. Defaulting to: '{output_filepath}'")
        else:
            output_to_console = True
            print("Output file could not be generated based on input extension. Results will be printed to console.")
    
    if not output_to_console: # Only print this if we're actually saving to a file
        print(f"Output file: '{output_filepath}'")


    # 3. Load the input table (DataFrame)
    print(f"Loading input data from: '{input_filepath}'...")
    try:
        ase_data = pd.read_csv(input_filepath, sep=separator)
        print(f"Successfully loaded {len(ase_data)} rows.")
    except Exception as e:
        print(f"Error loading input file '{input_filepath}': {e}", file=sys.stderr)
        sys.exit(1)

    # Check for required columns
    required_columns = [
        'gene_id', 'gene_name', 'gene_chr', 'gene_start', 'gene_end', 'gene_strand',
        'refCount', 'altCount', 'totalCount', 'contig', 'position', 'refAllele', 'altAllele']
    
    missing_columns = [col for col in required_columns if col not in ase_data.columns]
    
    if missing_columns:
        print(f"Error: Input file is missing the following required columns: {', '.join(missing_columns)}", file=sys.stderr)
        sys.exit(1)

    # 4. Call your gene_level_aggregation function
    try:
        gene_level_results = gene_level_aggregation(ase_data, threshold=threshold)
    except Exception as e:
        print(f"An error occurred during gene-level aggregation: {e}", file=sys.stderr)
        sys.exit(1)

    # 5. Handle output
    if not output_to_console and output_filepath:
        print(f"Saving results to: '{output_filepath}'")
        try:
            gene_level_results.to_csv(output_filepath, index=False, sep=separator) # Use determined separator for output
            print("Results saved successfully.")
        except Exception as e:
            print(f"Error writing results to output file '{output_filepath}': {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("\n--- Aggregated Gene-Level Results (Printed to Console) ---")
        print(gene_level_results.to_string()) # Use to_string() to print entire DataFrame
        print("\n--- End of Results ---")

    print(f"--- Gene Aggregation Complete ---")