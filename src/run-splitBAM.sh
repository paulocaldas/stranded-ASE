#!/bin/bash

# ---
# Script: run-splitBAM.sh
# Description: Splits a paired-end BAM file into forward and reverse strand-specific BAM files.
#              This is crucial for stranded RNA-seq data analysis (e.g., before ASEReadCounter).
#              The logic implemented here is typical for a "Reverse-Forward" (RF) stranded protocol
#              (e.g., Illumina TruSeq Stranded mRNA, NEBNext Ultra Directional).
#              Output BAM files are sorted and indexed.
# ---

# Function to display usage information
usage() {
  echo "Usage: $0 <input_bam> [options]"
  echo ""
  echo "Arguments:"
  echo "  <input_bam>      Path to the input BAM file (e.g., my_aligned.RG.bam)."
  echo "                   The output prefix will be inferred from this filename (e.g., 'my_aligned.RG')."
  echo ""
  echo "Options:"
  echo "  --threads <INT>  Number of threads for samtools (default: 1)."
  echo "  --output-dir <DIR> Directory to save final output BAMs (default: 'bams-split')."
  echo "  -h, --help       Display this help message."
  echo ""
  echo "Example:"
  echo "  $0 /path/to/my_sample.RG.bam --threads 8 --output-dir /data/split_bams"
  exit 1
}

# ---
# Parse command-line arguments
# ---

# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
fi

# Assign positional argument
INPUT_BAM="$1"
shift # Consume INPUT_BAM

# Validate required positional argument
if [[ -z "$INPUT_BAM" ]]; then
  echo "Error: Missing input BAM file."
  usage
fi

# Initialize optional parameters
THREADS=8 # Default threads
OUTPUT_DIR="bams-split" # Default output directory

# Parse optional arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --threads)
      THREADS="$2"
      shift # past argument
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift # past argument
      ;;
    *)
      echo "Error: Unknown option '$1'"
      usage
      ;;
  esac
  shift # past argument or value
done

# ---
# Validate input file existence
# ---
if [[ ! -f "$INPUT_BAM" ]]; then
  echo "Error: Input BAM file '${INPUT_BAM}' not found."
  exit 1
fi

# ---
# Infer output prefix from input BAM filename
# ---
# Example: "/path/to/my_sample.RG.bam" -> "my_sample.RG"
NAME_PREFIX=$(basename "${INPUT_BAM%.bam}")

# ---
# Create the output directory if it doesn't exist
# ---
echo "--- Starting BAM Splitting for ${INPUT_BAM} ---"
echo "  Inferred Output Prefix: ${NAME_PREFIX}"
echo "  Output Directory: ${OUTPUT_DIR}"
echo "  Threads: ${THREADS}"
echo "---"

mkdir -p "$OUTPUT_DIR" || { echo "Error: Could not create output directory '${OUTPUT_DIR}'."; exit 1; }
echo "  Ensured existence of folder: $OUTPUT_DIR"
echo "---"

# ---
# Forward strand logic (for RF stranded protocol):
# 1. Alignments of the second in pair if they map to the forward strand (0x80 set, 0x10 NOT set).
# 2. Alignments of the first in pair if their mate maps to the forward strand (0x40 set, 0x20 NOT set).
# ---

echo "## Processing Forward Strand Reads..."

echo "  Creating ${NAME_PREFIX}.fwd1.bam (Second in pair, mapped to forward strand)..."
# -f 128 (0x80): read is second in pair
# -F 16 (0x10): read is NOT mapped to reverse strand (i.e., mapped to forward strand)
samtools view -@ "$THREADS" -b -f 128 -F 16 "$INPUT_BAM" > "${NAME_PREFIX}.fwd1.bam" || \
  { echo "Error: Failed to create ${NAME_PREFIX}.fwd1.bam. Exiting."; exit 1; }
samtools index -@ "$THREADS" "${NAME_PREFIX}.fwd1.bam" || \
  { echo "Error: Failed to index ${NAME_PREFIX}.fwd1.bam. Exiting."; exit 1; }

echo "  Creating ${NAME_PREFIX}.fwd2.bam (First in pair, mate mapped to forward strand)..."
# -f 64 (0x40): read is first in pair
# -F 32 (0x20): mate is NOT mapped to reverse strand (i.e., mate mapped to forward strand)
samtools view -@ "$THREADS" -b -f 64 -F 32 "$INPUT_BAM" > "${NAME_PREFIX}.fwd2.bam" || \
  { echo "Error: Failed to create ${NAME_PREFIX}.fwd2.bam. Exiting."; exit 1; }
samtools index -@ "$THREADS" "${NAME_PREFIX}.fwd2.bam" || \
  { echo "Error: Failed to index ${NAME_PREFIX}.fwd2.bam. Exiting."; exit 1; }

# Combine alignments that originate on the forward strand.
echo "  Merging forward strand alignments into ${NAME_PREFIX}.fwd.bam ..."
samtools merge -@ "$THREADS" -f "${NAME_PREFIX}.fwd.bam" "${NAME_PREFIX}.fwd1.bam" "${NAME_PREFIX}.fwd2.bam" || \
  { echo "Error: Failed to merge forward strand BAMs. Exiting."; exit 1; }
samtools index -@ "$THREADS" "${NAME_PREFIX}.fwd.bam" || \
  { echo "Error: Failed to index ${NAME_PREFIX}.fwd.bam. Exiting."; exit 1; }
echo "  Forward strand BAM created: ${NAME_PREFIX}.fwd.bam"
echo "---"

# ---
# Reverse strand logic (for RF stranded protocol):
# 1. Alignments of the second in pair if they map to the reverse strand (0x80 set, 0x10 set).
# 2. Alignments of the first in pair if their mate maps to the reverse strand (0x40 set, 0x20 set).
# ---

echo "## Processing Reverse Strand Reads..."

echo "  Creating ${NAME_PREFIX}.rev1.bam (Second in pair, mapped to reverse strand)..."
# -f 144 (0x80 + 0x10): read is second in pair AND mapped to reverse strand
samtools view -@ "$THREADS" -b -f 144 "$INPUT_BAM" > "${NAME_PREFIX}.rev1.bam" || \
  { echo "Error: Failed to create ${NAME_PREFIX}.rev1.bam. Exiting."; exit 1; }
samtools index -@ "$THREADS" "${NAME_PREFIX}.rev1.bam" || \
  { echo "Error: Failed to index ${NAME_PREFIX}.rev1.bam. Exiting."; exit 1; }

echo "  Creating ${NAME_PREFIX}.rev2.bam (First in pair, mate mapped to reverse strand)..."
# -f 96 (0x40 + 0x20): read is first in pair AND mate is mapped to reverse strand
samtools view -@ "$THREADS" -b -f 96 "$INPUT_BAM" > "${NAME_PREFIX}.rev2.bam" || \
  { echo "Error: Failed to create ${NAME_PREFIX}.rev2.bam. Exiting."; exit 1; }
samtools index -@ "$THREADS" "${NAME_PREFIX}.rev2.bam" || \
  { echo "Error: Failed to index ${NAME_PREFIX}.rev2.bam. Exiting."; exit 1; }

# Combine alignments that originate on the reverse strand.
echo "  Merging reverse strand alignments into ${NAME_PREFIX}.rev.bam ..."
samtools merge -@ "$THREADS" -f "${NAME_PREFIX}.rev.bam" "${NAME_PREFIX}.rev1.bam" "${NAME_PREFIX}.rev2.bam" || \
  { echo "Error: Failed to merge reverse strand BAMs. Exiting."; exit 1; }
samtools index -@ "$THREADS" "${NAME_PREFIX}.rev.bam" || \
  { echo "Error: Failed to index ${NAME_PREFIX}.rev.bam. Exiting."; exit 1; }
echo "  Reverse strand BAM created: ${NAME_PREFIX}.rev.bam"
echo "---"

# ---
# Move final BAM files to the specified output directory
# ---
echo "## Moving final BAM files to $OUTPUT_DIR ..."
mv "${NAME_PREFIX}.fwd.bam" "$OUTPUT_DIR/" || \
  { echo "Error: Failed to move ${NAME_PREFIX}.fwd.bam. Exiting."; exit 1; }
mv "${NAME_PREFIX}.fwd.bam.bai" "$OUTPUT_DIR/" || \
  { echo "Error: Failed to move ${NAME_PREFIX}.fwd.bam.bai. Exiting."; exit 1; }
mv "${NAME_PREFIX}.rev.bam" "$OUTPUT_DIR/" || \
  { echo "Error: Failed to move ${NAME_PREFIX}.rev.bam. Exiting."; exit 1; }
mv "${NAME_PREFIX}.rev.bam.bai" "$OUTPUT_DIR/" || \
  { echo "Error: Failed to move ${NAME_PREFIX}.rev.bam.bai. Exiting."; exit 1; }
echo "---"

# ---
# Delete temporary files
# ---
echo "## Deleting temporary files ..."
rm "${NAME_PREFIX}.fwd1.bam" "${NAME_PREFIX}.fwd1.bam.bai" || \
  echo "Warning: Could not remove ${NAME_PREFIX}.fwd1.bam or its index."
rm "${NAME_PREFIX}.fwd2.bam" "${NAME_PREFIX}.fwd2.bam.bai" || \
  echo "Warning: Could not remove ${NAME_PREFIX}.fwd2.bam or its index."
rm "${NAME_PREFIX}.rev1.bam" "${NAME_PREFIX}.rev1.bam.bai" || \
  echo "Warning: Could not remove ${NAME_PREFIX}.rev1.bam or its index."
rm "${NAME_PREFIX}.rev2.bam" "${NAME_PREFIX}.rev2.bam.bai" || \
  echo "Warning: Could not remove ${NAME_PREFIX}.rev2.bam or its index."
echo "---"

echo "### Script completed successfully. ###"
