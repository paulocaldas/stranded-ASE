#!/bin/bash

# ---
# Script for annotating ASEReadCounter .tab files with gene information from a GTF.
# Steps:
# 0. Accept two .tab files (.fwd.tab, .rev.tab) and a GTF file as input.
# 1. Create a 'tmp' directory for intermediate files.
# 2. Convert GTF to BED format and save in 'tmp'.
# 3. Split the BED GTF by strand ('+' and '-') into 'tmp/fwd.bed' and 'tmp/rev.bed'.
# 4. Convert input .tab files to temporary BED format for bedtools intersect.
# 5. Use bedtools intersect to cross each input tab (now temp BED) with the corresponding GTF BED file.
# 6. Reconstruct the output with headers and save in the input tables' directory.
# ---

# Function to display usage information
usage() {
  echo "Usage: $0 <fwd_tab_file> <rev_tab_file> <gtf_file>"
  echo ""
  echo "Arguments:"
  echo "  <fwd_tab_file>   Path to the ASEReadCounter forward strand .tab file (e.g., my_sample.fwd.tab)"
  echo "  <rev_tab_file>   Path to the ASEReadCounter reverse strand .tab file (e.g., my_sample.rev.tab)"
  echo "  <gtf_file>       Path to the GTF annotation file (e.g., genes.gtf)"
  echo ""
  echo "Example:"
  echo "  $0 ase_output/sample1.fwd.tab ase_output/sample1.rev.tab annotations/genes.gtf"
  exit 1
}

# ---
# Parse command-line arguments
# ---

# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
fi

# Assign positional arguments
FWD_TAB_FILE="$1"
REV_TAB_FILE="$2"
GTF_FILE="$3"

# Validate required arguments
if [[ -z "$FWD_TAB_FILE" || -z "$REV_TAB_FILE" || -z "$GTF_FILE" ]]; then
  echo "Error: Missing required arguments."
  usage
fi

# Validate input files exist
if [[ ! -f "$FWD_TAB_FILE" ]]; then echo "Error: Forward tab file '${FWD_TAB_FILE}' not found."; exit 1; fi
if [[ ! -f "$REV_TAB_FILE" ]]; then echo "Error: Reverse tab file '${REV_TAB_FILE}' not found."; exit 1; fi
if [[ ! -f "$GTF_FILE" ]]; then echo "Error: GTF file '${GTF_FILE}' not found."; exit 1; fi

# Define paths and variables
TAB_DIR=$(dirname "${FWD_TAB_FILE}")
FWD_BASE_NAME=$(basename "${FWD_TAB_FILE%.tab}") # e.g., "my_sample.fwd"
REV_BASE_NAME=$(basename "${REV_TAB_FILE%.tab}") # e.g., "my_sample.rev"

# GTF file base name
GTF_BASE_NAME=$(basename "${GTF_FILE%.gtf}")

# Temporary directory
TMP_DIR="${TAB_DIR}/tmp" # Place tmp folder within the same dir as the input tabs

# ---
# Pipeline Steps
# ---

echo "--- Starting Annotation Pipeline ---"
echo "  Forward TAB: ${FWD_TAB_FILE}"
echo "  Reverse TAB: ${REV_TAB_FILE}"
echo "  GTF File: ${GTF_FILE}"
echo "  Temporary directory: ${TMP_DIR}"
echo "---"

# Step 1: Create 'tmp' directory
echo ".creating temporary directory"
if [ ! -d "$TMP_DIR" ]; then
  mkdir -p "$TMP_DIR" || { echo "Error: Could not create temporary directory."; exit 1; }
  echo "  Directory '${TMP_DIR}' created."
else
  echo "  Directory '${TMP_DIR}' already exists. Skipping creation."
fi
echo " "

# Step 2: Convert GTF (gene features) to BED format
echo ".converting GTF to BED format"
GTF_BED="${TMP_DIR}/${GTF_BASE_NAME}.bed"
if [ ! -f "$GTF_BED" ]; then
  cat "$GTF_FILE" | \
    awk 'OFS="\t" {
      if ($3=="gene") {
        gsub("gene_id ", "", $10);  # Clean gene_id
        gsub("gene_name ", "", $14); # Clean gene_name
        print $1, $4-1, $5, $10, $14, $7
      }
    }' | tr -d '";' > "${GTF_BED}" || \
    { echo "Error: GTF to BED conversion failed."; exit 1; }
  echo "  GTF converted to ${GTF_BED}."
else
  echo "  GTF BED file '${GTF_BED}' already exists. Skipping conversion."
fi
echo " "

# Step 3: Split GTF BED by strand
echo ".splitting GTF BED by strand"
GTF_FWD_BED="${TMP_DIR}/${GTF_BASE_NAME}.fwd.bed"
GTF_REV_BED="${TMP_DIR}/${GTF_BASE_NAME}.rev.bed"

if [ ! -f "$GTF_FWD_BED" ] || [ ! -f "$GTF_REV_BED" ]; then
  awk 'OFS="\t" {if ($6=="+") print}' "${GTF_BED}" > "${GTF_FWD_BED}" || \
    { echo "Error: Splitting GTF BED (forward) failed."; exit 1; }
  awk 'OFS="\t" {if ($6=="-") print}' "${GTF_BED}" > "${GTF_REV_BED}" || \
    { echo "Error: Splitting GTF BED (reverse) failed."; exit 1; }
  echo "  GTF BED split into ${GTF_FWD_BED} and ${GTF_REV_BED}."
else
  echo "  Strand-specific GTF BED files already exist. Skipping splitting."
fi
echo " "

# Step 4: Convert Input .tab files to "rich" temporary BED format for intersect
echo ".converting input .tab files to temporary BED format (preserving all columns)"
TEMP_FWD_BED="${TMP_DIR}/${FWD_BASE_NAME}.temp.bed"
TEMP_REV_BED="${TMP_DIR}/${REV_BASE_NAME}.temp.bed"

awk 'NR > 1 {
    # Store original columns to append later
    orig_cols = $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9

    # BED core fields
    chrom = $1
    chromStart = ($2 - 1)
    chromEnd = $2
    name = ($3 == "." ? $1 "_" $2 "_" $4 "_" $5 : $3)
    score = $8 # totalCount
    strand = "."

    print chrom "\t" chromStart "\t" chromEnd "\t" name "\t" score "\t" strand "\t" orig_cols
}' "${FWD_TAB_FILE}" > "${TEMP_FWD_BED}" || \
{ echo "Error: Converting FWD .tab to temp BED failed."; exit 1; }
echo "  Created temporary BED for forward reads: ${TEMP_FWD_BED}"

awk 'NR > 1 {
    orig_cols = $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9

    chrom = $1
    chromStart = ($2 - 1)
    chromEnd = $2
    name = ($3 == "." ? $1 "_" $2 "_" $4 "_" $5 : $3)
    score = $8 # totalCount
    strand = "."

    print chrom "\t" chromStart "\t" chromEnd "\t" name "\t" score "\t" strand "\t" orig_cols
}' "${REV_TAB_FILE}" > "${TEMP_REV_BED}" || \
{ echo "Error: Converting REV .tab to temp BED failed."; exit 1; }
echo "  Created temporary BED for reverse reads: ${TEMP_REV_BED}"
echo " "

# Step 5: Perform bedtools intersect
echo ".performing bedtools intersect"

# Final output files will be in the same directory as input tabs
ANNOTATED_FWD_FILE="${TAB_DIR}/${FWD_BASE_NAME}.annotated.tab"
ANNOTATED_REV_FILE="${TAB_DIR}/${REV_BASE_NAME}.annotated.tab"

# Construct the new header for the annotated output
NEW_HEADER="contig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\ttotalCount\tsample\tgene_chr\tgene_start\tgene_end\tgene_id\tgene_name\tgene_strand"

echo -e "${NEW_HEADER}" > "${ANNOTATED_FWD_FILE}" # Write header to fwd output
echo -e "${NEW_HEADER}" > "${ANNOTATED_REV_FILE}" # Write header to rev output

echo "  Intersecting forward reads with forward GTF genes..."
bedtools intersect -a "${TEMP_FWD_BED}" -b "${GTF_FWD_BED}" -wa -wb | \
  awk 'OFS="\t" {
    print $1, ($2 + 1), $7, $8, $9, $10, $11, $12, $13, \
          $14, $15, $16, $17, $18, $19
  }' >> "${ANNOTATED_FWD_FILE}" || \
  { echo "Error: bedtools intersect for forward reads failed."; exit 1; }
echo "  Annotated forward reads saved to: ${ANNOTATED_FWD_FILE}"

echo "  Intersecting reverse reads with reverse GTF genes..."
bedtools intersect -a "${TEMP_REV_BED}" -b "${GTF_REV_BED}" -wa -wb | \
  awk 'OFS="\t" {
    print $1, ($2 + 1), $7, $8, $9, $10, $11, $12, $13, \
          $14, $15, $16, $17, $18, $19
  }' >> "${ANNOTATED_REV_FILE}" || \
  { echo "Error: bedtools intersect for reverse reads failed."; exit 1; }
echo "  Annotated reverse reads saved to: ${ANNOTATED_REV_FILE}"
echo " "

# Step 6: Merge FWD and REV files into a single table
echo ".merging fwd and rev tables"

OUTPUT_MERGED_FILE=${ANNOTATED_FWD_FILE/.fwd/}
echo "Output File: $OUTPUT_MERGED_FILE"

# Get the header from the first file and write it to the output file
head -n 1 "$ANNOTATED_REV_FILE" > "$OUTPUT_MERGED_FILE" || \
  { echo "Error: Failed to write header from ${FWD_ANNOTATED_FILE}. Exiting."; exit 1; }

# Append the data (excluding header) from the first file
tail -n +2 "$ANNOTATED_FWD_FILE" >> "$OUTPUT_MERGED_FILE" || \
  { echo "Error: Failed to append data from ${FWD_ANNOTATED_FILE}. Exiting."; exit 1; }

# 3. Append the data (excluding header) from the second file
tail -n +2 "$ANNOTATED_REV_FILE" >> "$OUTPUT_MERGED_FILE" || \
  { echo "Error: Failed to append data from ${REV_ANNOTATED_FILE}. Exiting."; exit 1; }

# Step 7: Clean up temporary files
echo ".cleaning up temporary files"
rm -rf "$TMP_DIR"
rm "$ANNOTATED_REV_FILE"
rm "$ANNOTATED_FWD_FILE"
echo " "

echo " "
echo "> Variant Annotation Completed!"
