# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 13:51:41 2025

@author: paulo
"""

import pandas as pd
import numpy as np
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

def gene_level_aggregation(ase_res, threshold=10):

    # Filter for minimum read number
    ase_res_filtered = ase_res[ase_res['totalCount'] > threshold].copy()

    if ase_res_filtered.empty:
        print(f"Warning: No SNPs left after applying threshold > {threshold}.")
        return None

    # Create variant label
    variant_label = (ase_res_filtered['contig'] + '_' +
                     ase_res_filtered['position'].astype(str) + '_' +
                     ase_res_filtered['refAllele'] + '_' +
                     ase_res_filtered['altAllele'])

    ase_res_filtered.insert(ase_res_filtered.shape[1], 'variant', variant_label)

    # Compute minor allele frequency for each SNP
    ase_res_filtered['maf'] = ase_res_filtered[['refCount', 'altCount']].min(axis=1) / ase_res_filtered['totalCount']

    # Weighted average of individual SNP MAFs
    def weighted_aggregation(df):
        # Drop groupby keys to avoid include_groups deprecation (pandas >= 2.0)
        df = df.drop(columns=['gene_id', 'gene_name', 'gene_chr',
                               'gene_start', 'gene_end', 'gene_strand'],
                     errors='ignore')
        ref_sum = df['refCount'].sum()
        alt_sum = df['altCount'].sum()
        total_sum = df['totalCount'].sum()
        weights = df['totalCount']
        maf = df['maf']
        weighted_maf = (maf * weights).sum() / weights.sum() if weights.sum() > 0 else 0
        variants_concat = ','.join(df['variant'].astype(str))

        return pd.Series({
            'refCount': ref_sum,
            'altCount': alt_sum,
            'totalCount': total_sum,
            'maf_weighted': weighted_maf,
            'variant': variants_concat})

    print("Aggregating SNP data to gene level with weighted MAF...")

    gene_level = (
        ase_res_filtered.groupby(['gene_id', 'gene_name', 'gene_chr', 'gene_start', 'gene_end', 'gene_strand'])
        .apply(weighted_aggregation)
        .reset_index())

    # Compute gene-level MAF based on summed counts
    gene_level.insert(gene_level.shape[1]-1, 'maf_tot',
                      gene_level[['refCount', 'altCount']].min(axis=1) / gene_level['totalCount'])

    # Add allelic expression label
    gene_level.insert(gene_level.shape[1]-1, 'allelic_exp',
                      gene_level['maf_tot'].apply(isBiallelic))

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
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'input_file',
        type=str,
        help="Path to the input file (CSV or TSV) containing SNP-level ASE results."
             "\nExpected columns:\n"
             "  - 'gene_id', 'gene_name', 'gene_chr', 'gene_start', 'gene_end', 'gene_strand'\n"
             "  - 'refCount', 'altCount', 'totalCount'\n"
             "  - 'contig', 'position', 'refAllele', 'altAllele'"
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help="Optional: Path to the output CSV/TSV file. "
             "If not provided, defaults to '<input_filename>.geneLevelMAF.<ext>'. "
             "If extension is unrecognized, results are printed to console."
    )

    parser.add_argument(
        '-t', '--threshold',
        type=int,
        default=10,
        help="Optional: Minimum total read count for SNPs to be included. "
             "SNPs with 'totalCount' <= threshold will be excluded. Default: 10."
    )

    args = parser.parse_args()

    input_filepath = args.input_file
    output_filepath = args.output
    threshold = args.threshold

    print(f"--- Starting Gene Aggregation ---")
    print(f"Input file: '{input_filepath}'")
    print(f"Threshold for totalCount: {threshold}")

    # 1. Validate input file
    if not os.path.exists(input_filepath):
        print(f"Error: Input file not found at '{input_filepath}'", file=sys.stderr)
        sys.exit(1)

    # Determine separator from file extension
    input_dir = os.path.dirname(input_filepath)
    input_basename = os.path.basename(input_filepath)
    filename_without_ext, file_extension = os.path.splitext(input_basename)
    file_extension = file_extension.lower()

    if file_extension == '.csv':
        separator = ','
    elif file_extension in ('.tab', '.tsv', '.txt'):
        separator = '\t'
    else:
        print(f"Warning: Input file format '{file_extension}' not explicitly supported. "
              "Attempting to read as CSV.", file=sys.stderr)
        separator = ','

    # 2. Determine output file path
    output_to_console = False
    if output_filepath is None:
        if file_extension in ('.tab', '.csv', '.tsv', '.txt'):
            default_output_filename = f"{filename_without_ext}.geneLevelMAF{file_extension}"
            output_filepath = os.path.join(input_dir, default_output_filename)
            print(f"No output file specified. Defaulting to: '{output_filepath}'")
        else:
            output_to_console = True
            print("Output file could not be determined. Results will be printed to console.")

    if not output_to_console:
        print(f"Output file: '{output_filepath}'")

    # 3. Load input table
    print(f"Loading input data from: '{input_filepath}'...")
    try:
        ase_data = pd.read_csv(input_filepath, sep=separator)
        print(f"Successfully loaded {len(ase_data)} rows.")
    except Exception as e:
        print(f"Error loading input file '{input_filepath}': {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Validate required columns
    required_columns = [
        'gene_id', 'gene_name', 'gene_chr', 'gene_start', 'gene_end', 'gene_strand',
        'refCount', 'altCount', 'totalCount', 'contig', 'position', 'refAllele', 'altAllele']

    missing_columns = [col for col in required_columns if col not in ase_data.columns]
    if missing_columns:
        print(f"Error: Input file is missing required columns: {', '.join(missing_columns)}", file=sys.stderr)
        sys.exit(1)

    # 5. Run gene-level aggregation
    try:
        gene_level_results = gene_level_aggregation(ase_data, threshold=threshold)
    except Exception as e:
        print(f"An error occurred during gene-level aggregation: {e}", file=sys.stderr)
        sys.exit(1)

    # 6. Check for empty results
    if gene_level_results is None:
        print("Error: No results to save. All SNPs were filtered out by the threshold.", file=sys.stderr)
        sys.exit(1)

    # 7. Handle output
    if not output_to_console and output_filepath:
        print(f"Saving results to: '{output_filepath}'")
        try:
            gene_level_results.to_csv(output_filepath, index=False, sep=separator)
            print("Results saved successfully.")
        except Exception as e:
            print(f"Error writing results to '{output_filepath}': {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("\n--- Aggregated Gene-Level Results ---")
        print(gene_level_results.to_string())
        print("\n--- End of Results ---")

    print(f"--- Gene Aggregation Complete ---")