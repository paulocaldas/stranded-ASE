#!/bin/bash

# ---
# Pipeline script to process a single BAM file for Allele-Specific Expression (ASE) analysis.
# This script performs:
# 0. Adding/replacing Read Groups (using run-addReadGroups.sh) -> NOW A PREPROCESSING STEP 
# 1. Splitting BAM by strand (for RNA-seq, using run-splitBams.sh)
# 2. Running GATK ASEReadCounter on both forward and reverse strands (using run-aseReadCounter.sh)
# 3. Annotate ASEReadCounter Output file (using run-annotateVariants.sh)
# 4. Compute gene-level Minor Allele Frequency (MAF)
# ---

# Function to display usage information
usage() {
  echo "Usage: $0 <input_bam> <variants_vcf> <reference_fasta> <referece_gtf> [options]"
  echo ""
  echo "Arguments:"
  echo "  <input_bam>        Path to the input BAM file containing proper read groups(e.g., /path/to/raw_reads.bam)"
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
  echo "  $0 data/my_sample.bam variants.vcf.gz ref.fasta ref.gtf --min-mapping-quality 30"
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
VARIANTS_VCF="$2"
REFERENCE_FASTA="$3"
REFERENCE_GTF="$4"

# Shift arguments to parse optional ASEReadCounter flags
shift 4

# Validate required positional arguments
if [[ -z "$INPUT_BAM" || -z "$VARIANTS_VCF" || -z "$REFERENCE_FASTA" || -z "$REFERENCE_GTF" ]]; then
  echo " "
  echo "Error: Missing required arguments."
  usage
fi

# Check if the input BAM file has Read Groups (@RG) in its header
echo ""
echo "Checking for Read Groups in ${INPUT_BAM}..."
if ! samtools view -H "${INPUT_BAM}" | grep -q '^@RG'; then
  echo " "
  echo "ERROR: Input BAM file '${INPUT_BAM}' does not contain any Read Group (@RG) information in its header." >&2
  echo "Read groups are required for downstream analysis (e.g., ASEReadCounter)." >&2
  echo "Please add read groups to your BAM file using Picard's AddOrReplaceReadGroups or a similar tool." >&2
  echo " "
  exit 1
fi
echo "Read Groups found in ${INPUT_BAM}."


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
#RG_BAM="${INPUT_BAM_DIR}/${BASE_FILENAME}.RG.bam" # outdated - when addReadGroups was part of the pipeline
RG_BAM="$INPUT_BAM"

# Output directory from run-splitBams.sh
# Assumes run-splitBams.sh creates a subdirectory named "bams-split"
SPLIT_BAM_SUBDIR="bams-split/"
FWD_BAM="${SPLIT_BAM_SUBDIR}/${BASE_FILENAME}.fwd.bam"
REV_BAM="${SPLIT_BAM_SUBDIR}/${BASE_FILENAME}.rev.bam"

# Output names for ASEReadCounter
ASE_COUNTER_SUBDIR="out-aseReadCounter/"
FWD_OUTPUT_NAME="${ASE_COUNTER_SUBDIR}/${BASE_FILENAME}.fwd.tab"
REV_OUTPUT_NAME="${ASE_COUNTER_SUBDIR}/${BASE_FILENAME}.rev.tab"

# output names for run-annotateVariants.sh
# already defined from above

# ---
# Pipeline Steps
# ---

echo "=== Starting ASE Pipeline for ${INPUT_BAM} ==="
echo "- VCF: ${VARIANTS_VCF}"
echo "- Reference: ${REFERENCE_FASTA}"
echo "- ASEReadCounter options: ${ASER_OPTIONS:-None}" # Show options if provided
echo "==="

# -- OLD --- Step 0: Add Read Groups to BAM file - NOW A PREPROCESSING STEP
#echo "# Step 1: Adding/Replacing Read Groups"
#echo "- Input BAM: ${INPUT_BAM}"
#echo "- VCF SAMPLE: ${VCF_SAMPLE_NAME}"
#echo "- Output: ${RG_BAM}"

#bash $SRC_CODE_DIR/run-addReadGroups.sh "${INPUT_BAM}" "${VCF_SAMPLE_NAME}" || \
#  { echo "Error: Step 1 (Add Read Groups) failed. Exiting."; exit 1; }

#echo ".Read groups added successfully."
#echo " "

# Step 1: Split BAM file by strand
echo "# Step 1: Splitting BAM by Strand"
echo "- Input: ${RG_BAM}"
echo "- Output Directory: ${SPLIT_BAM_SUBDIR}"
bash $SRC_CODE_DIR/run-splitBAM.sh "${RG_BAM}" || \
  { echo "Error: Step 2 (Split BAM) failed. Exiting."; exit 1; }

echo ".BAM split successfully to ${SPLIT_BAM_SUBDIR}."
echo " "

# Step 2A: Run ASEReadCounter on Forward Reads
echo "# Step 2A: Running ASEReadCounter on Forward Reads"
echo "- Input: ${FWD_BAM}"
echo "- Output Name: ${FWD_OUTPUT_NAME}"
bash $SRC_CODE_DIR/run-aseReadCounter.sh "${FWD_BAM}" "${VARIANTS_VCF}" "${REFERENCE_FASTA}" \
  --output-name "${FWD_OUTPUT_NAME}" ${ASER_OPTIONS} || \
  { echo "Error: Step 3a (ASEReadCounter - Forward) failed. Exiting."; exit 1; }

echo ".ASEReadCounter completed for forward reads."
echo " "

# Step 2B: Run ASEReadCounter on Reverse Reads
echo "# Step 2B: Running ASEReadCounter on Reverse Reads"
echo "- Input: ${REV_BAM}"
echo "- Output Name: ${REV_OUTPUT_NAME}"

bash $SRC_CODE_DIR/run-aseReadCounter.sh "${REV_BAM}" "${VARIANTS_VCF}" "${REFERENCE_FASTA}" \
  --output-name "${REV_OUTPUT_NAME}" ${ASER_OPTIONS} || \
  { echo "Error: Step 3b (ASEReadCounter - Reverse) failed. Exiting."; exit 1; }

echo ".ASEReadCounter completed for reverse reads."

echo ".Delete Intermediate Files in ${SPLIT_BAM_SUBDIR}"
rm -r "${SPLIT_BAM_SUBDIR}"
echo " "

# Step 3: Annotate Output from ASEReadCounter
echo "# Step 3: Annotate ASEReadCounter Output"
echo "- Input 1: $FWD_OUTPUT_NAME"
echo "- Input 2: $REV_OUTPUT_NAME"
echo "- GTF File: $REFERENCE_GTF"

bash $SRC_CODE_DIR/run-annotateVariants.sh "$FWD_OUTPUT_NAME" "$REV_OUTPUT_NAME" "$REFERENCE_GTF"

# Step 4: Compute MAF for each gene
OUTPUT_MERGED_FILE="${ASE_COUNTER_SUBDIR}/${BASE_FILENAME}.annotated.tab"

echo " Step 4: Compute MAF for each Gene"
echo "- Input: $OUTPUT_MERGED_FILE"

python $SRC_CODE_DIR/run-geneLevelMAF.py $OUTPUT_MERGED_FILE

# The End
echo "===  Pipeline completed successfully for ${INPUT_BAM} ==="
