#!/bin/bash

# ---
# Pipeline script to process a single BAM file for Allele-Specific Expression (ASE) analysis.
# This script performs:
# 1. Adding/replacing Read Groups (using run-addReadGroups.sh)
# 2. Splitting BAM by strand (for RNA-seq, using run-splitBams.sh)
# 3. Running GATK ASEReadCounter on both forward and reverse strands (using run-aseReadCounter.sh)
# 4. Annotate ASEReadCounter Output file (using run-annotateVariants.sh)
# ---

# Function to display usage information
usage() {
  echo "Usage: $0 <input_bam> <vcf_sample_name> <variants_vcf> <reference_fasta> [options]"
  echo ""
  echo "Arguments:"
  echo "  <input_bam>        Path to the input BAM file (e.g., /path/to/raw_reads.bam)"
  echo "  <vcf_sample_name>  The biological sample name (e.g., 'hESC_H9'). This MUST match a sample column in your VCF."
  echo "  <variants_vcf>     Path to the input VCF file with variant sites (e.g., variants.vcf.gz)"
  echo "  <reference_fasta>  Path to the reference FASTA file (e.g., genome.fa)"
  echo "  <reference_gtf>    Path to the reference GTF file (e.g., genome.gtf)"
  echo ""
  echo "Options (for ASEReadCounter, pass-through to run-aseReadCounter.sh):"
  echo "  --min-mapping-quality <INT> Minimum mapping quality (MQ) for reads (default: 20)"
  echo "  --min-base-quality <INT>  Minimum base quality (BQ) for bases (default: 20)"
  echo "  --perform-indels   Include indels in the counting (flag)"
  echo "  --verbosity <LEVEL> Set GATK verbosity level (default: INFO)"
  echo "  -h, --help         Display this help message"
  echo ""
  echo "Example:"
  echo "  $0 data/my_sample.bam hESC_H9 variants.vcf.gz ref.fasta ref.gtf --min-mapping-quality 30"
  exit 1
}

# ---
# Parse command-line arguments
# ---

# Check for help flag first
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
fi

# Assign positional arguments
INPUT_BAM="$1"
VCF_SAMPLE_NAME="$2"
VARIANTS_VCF="$3"
REFERENCE_FASTA="$4"
REFERENCE_GTF="$5"

# Shift arguments to parse optional ASEReadCounter flags
shift 5

# Validate required positional arguments
if [[ -z "$INPUT_BAM" || -z "$VCF_SAMPLE_NAME" || -z "$VARIANTS_VCF" || -z "$REFERENCE_FASTA" || -z "$REFERENCE_GTF" ]]; then
  echo " "
  echo "Error: Missing required arguments."
  usage
fi

# Derive base filename and directory for intermediate files
# This ensures paths are correctly handled regardless of input_bam's location
SRC_CODE_DIR="src" # source code directory
BASE_FILENAME=$(basename "${INPUT_BAM%.bam}") # e.g., "my_sample" from "/path/to/my_sample.bam"
INPUT_BAM_DIR=$(dirname "$INPUT_BAM")         # e.g., "/path/to"

# Collect any additional options to pass to ASEReadCounter
# This captures all remaining arguments like --min-mapping-quality
ASER_OPTIONS="$@"

# ---
# Define Paths for Intermediate and Output Files
# ---

# Output from run-addReadGroups.sh
RG_BAM="${INPUT_BAM_DIR}/${BASE_FILENAME}.RG.bam"

# Output directory from run-splitBams.sh
# Assumes run-splitBams.sh creates a subdirectory named "bams-split"
SPLIT_BAM_SUBDIR="bams-split/"
FWD_BAM="${SPLIT_BAM_SUBDIR}/${BASE_FILENAME}.RG.fwd.bam"
REV_BAM="${SPLIT_BAM_SUBDIR}/${BASE_FILENAME}.RG.rev.bam"

# Output names for ASEReadCounter
ASE_COUNTER_SUBDIR="out-aseReadCounter/"
FWD_OUTPUT_NAME="${ASE_COUNTER_SUBDIR}/${BASE_FILENAME}.fwd.tab"
REV_OUTPUT_NAME="${ASE_COUNTER_SUBDIR}/${BASE_FILENAME}.rev.tab"

# output names for run-annotateVariants.sh
# already defined from above

# ---
# Pipeline Steps
# ---

echo "--- Starting ASE Pipeline for ${INPUT_BAM} (VCF Sample: ${VCF_SAMPLE_NAME}) ---"
echo " VCF: ${VARIANTS_VCF}"
echo " Reference: ${REFERENCE_FASTA}"
echo " ASEReadCounter options: ${ASER_OPTIONS:-None}" # Show options if provided
echo "---"

# Step 1: Add Read Groups to BAM file
echo "# Step 1: Adding/Replacing Read Groups"
echo " Input BAM: ${INPUT_BAM}"
echo " VCF SAMPLE: ${VCF_SAMPLE_NAME}"
echo " Output: ${RG_BAM}"

bash $SRC_CODE_DIR/run-addReadGroups.sh "${INPUT_BAM}" "${VCF_SAMPLE_NAME}" || \
  { echo "Error: Step 1 (Add Read Groups) failed. Exiting."; exit 1; }

echo "Read groups added successfully."
echo " "

# Step 2: Split BAM file by strand
echo "# Step 2: Splitting BAM by Strand"
echo "  Input: ${RG_BAM}"
echo "  Output Directory: ${SPLIT_BAM_SUBDIR}"
bash $SRC_CODE_DIR/run-splitBAM.sh "${RG_BAM}" || \
  { echo "Error: Step 2 (Split BAM) failed. Exiting."; exit 1; }

echo "BAM split successfully."

# delete duplicated bam from initial folder
rm ${RG_BAM}
rm ${RG_BAM%.*}.bai
echo "Duplicated BAMs - containing ReadGroups - where deleted"
echo " "

# Step 3A: Run ASEReadCounter on Forward Reads
echo "# Step 3A: Running ASEReadCounter on Forward Reads"
echo "  Input: ${FWD_BAM}"
echo "  Output Name: ${FWD_OUTPUT_NAME}"
bash $SRC_CODE_DIR/run-aseReadCounter.sh "${FWD_BAM}" "${VARIANTS_VCF}" "${REFERENCE_FASTA}" \
  --output-name "${FWD_OUTPUT_NAME}" ${ASER_OPTIONS} || \
  { echo "Error: Step 3a (ASEReadCounter - Forward) failed. Exiting."; exit 1; }

echo "ASEReadCounter completed for forward reads."
echo " "

# Step 3A: Run ASEReadCounter on Reverse Reads
echo "## Step 3B: Running ASEReadCounter on Reverse Reads"
echo "  Input: ${REV_BAM}"
echo "  Output Name: ${REV_OUTPUT_NAME}"

bash $SRC_CODE_DIR/run-aseReadCounter.sh "${REV_BAM}" "${VARIANTS_VCF}" "${REFERENCE_FASTA}" \
  --output-name "${REV_OUTPUT_NAME}" ${ASER_OPTIONS} || \
  { echo "Error: Step 3b (ASEReadCounter - Reverse) failed. Exiting."; exit 1; }

echo "ASEReadCounter completed for reverse reads."
echo " "

# Step 4: Annotate Output from ASEReadCounter
echo "Step 4: Annotate ASEReadCounter Output"
echo "Input 1: $FWD_OUTPUT_NAME"
echo "Input 2: $REV_OUTPUT_NAME"
echo "GTF File: $REFERENCE_GTF"

bash $SRC_CODE_DIR/run-annotateVariants.sh "$FWD_OUTPUT_NAME" "$REV_OUTPUT_NAME" "$REFERENCE_GTF"

# Combine FWD and REV annotation into one file and remove intermediate files

# The End
echo "# Pipeline completed successfully for ${INPUT_BAM} #"
