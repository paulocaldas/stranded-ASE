#!/bin/bash

# ---
# Script to run GATK ASEReadCounter on a single BAM file.
# This script expects to be called with a single BAM input (e.g., a strand-specific BAM)
# and required VCF, Reference FASTA files, and sample name.
# ---

# Function to display usage information
usage() {
  echo "Usage: $0 <input_bam> <variants_vcf> <reference_fasta> <sample_name> [options]"
  echo ""
  echo "Arguments:"
  echo "  <input_bam>        Path to the input BAM file (e.g., my_sample.fwd.bam or my_sample.rev.bam)"
  echo "  <variants_vcf>     Path to the input VCF file with variant sites (e.g., heterozygous_snps.vcf.gz)"
  echo "  <reference_fasta>  Path to the reference FASTA file (e.g., genome.fa)"
  echo "  <sample_name>      Sample name as in VCF genotype columns (e.g., hESC_H9)"
  echo ""
  echo "Options:"
  echo "  --output-name <file>  Name for the output .tab file (default: <input_basename>.tab in 'out-aseReadCounter' folder)"
  echo "  --min-mapping-quality <INT> Minimum mapping quality (MQ) for reads (default: 20)"
  echo "  --min-base-quality <INT>  Minimum base quality (BQ) for bases (default: 20)"
  echo "  --perform-indels   Include indels in the counting (flag, use if VCF includes indels)"
  echo "  --verbosity <LEVEL> Set GATK verbosity level (default: INFO, options: ERROR, WARNING, INFO, DEBUG, TRACE)"
  echo "  -h, --help            Display this help message"
  echo ""
  echo "Example:"
  echo "  $0 bams/SRR123.fwd.bam variants/SRR123.het.vcf.gz ref.fasta hESC_H9 --output-name SRR123.fwd.ase.tab --min-mapping-quality 30"
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
INPUT_BAM="$1"
VARIANTS_VCF="$2"
REFERENCE_FASTA="$3"
SAMPLE_NAME="$4"

# Shift positional arguments
shift 4

# Validate required positional arguments
if [[ -z "$INPUT_BAM" || -z "$VARIANTS_VCF" || -z "$REFERENCE_FASTA" || -z "$SAMPLE_NAME" ]]; then
  echo "Error: Missing required arguments."
  usage
fi

# Initialize optional parameters
OUTPUT_NAME_ARG=""
MIN_MQ="20"
MIN_BQ="20"
PERFORM_INDELS_FLAG=""
VERBOSITY_LEVEL="INFO"

# Parse optional arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --output-name)
      OUTPUT_NAME_ARG="$2"
      shift
      ;;
    --min-mapping-quality)
      MIN_MQ="$2"
      shift
      ;;
    --min-base-quality)
      MIN_BQ="$2"
      shift
      ;;
    --perform-indels)
      PERFORM_INDELS_FLAG="--perform-indels"
      ;;
    --verbosity)
      VERBOSITY_LEVEL="$2"
      shift
      ;;
    *)
      echo "Error: Unknown option '$1'"
      usage
      ;;
  esac
  shift
done

# ---
# Validate input files
# ---
if [[ ! -f "$INPUT_BAM" ]]; then echo "Error: Input BAM file '${INPUT_BAM}' not found."; exit 1; fi
if [[ ! -f "$VARIANTS_VCF" ]]; then echo "Error: Variants VCF file '${VARIANTS_VCF}' not found."; exit 1; fi
if [[ ! -f "$REFERENCE_FASTA" ]]; then echo "Error: Reference FASTA file '${REFERENCE_FASTA}' not found."; exit 1; fi

# ---
# Determine final output file path
# ---
FINAL_OUTPUT_NAME=""
DEFAULT_OUTPUT_DIR="out-aseReadCounter"

if [[ -z "$OUTPUT_NAME_ARG" ]]; then
  mkdir -p "$DEFAULT_OUTPUT_DIR" || { echo "Error: Could not create default output directory '${DEFAULT_OUTPUT_DIR}'."; exit 1; }
  BASE_FILENAME=$(basename "${INPUT_BAM%.bam}")
  FINAL_OUTPUT_NAME="${DEFAULT_OUTPUT_DIR}/${BASE_FILENAME}.tab"
else
  FINAL_OUTPUT_NAME="${OUTPUT_NAME_ARG}"
  OUTPUT_DIR_PATH=$(dirname "${FINAL_OUTPUT_NAME}")
  if [[ -n "$OUTPUT_DIR_PATH" && "$OUTPUT_DIR_PATH" != "." ]]; then
    mkdir -p "$OUTPUT_DIR_PATH" || { echo "Error: Could not create output directory for '${FINAL_OUTPUT_NAME}'."; exit 1; }
  fi
fi

# ---
# Filter VCF file to contain only Heterozygous sites for SAMPLE_NAME
# ---
VARIANTS_HET_VCF="${VARIANTS_VCF/.vcf.gz/.HET.vcf.gz}"

if [ ! -f "${VARIANTS_HET_VCF}" ]; then
  echo ""
  echo "Creating heterozygous-only VCF for sample '${SAMPLE_NAME}': ${VARIANTS_HET_VCF}"
  echo ""

  gatk SelectVariants \
    -V "${VARIANTS_VCF}" \
    --select "vc.getGenotype('${SAMPLE_NAME}').isHet()" \
    -O "${VARIANTS_HET_VCF}" || { echo "Error: SelectVariants failed."; exit 1; }
else
  echo ""
  echo "Heterozygous VCF already exists: ${VARIANTS_HET_VCF}, skipping."
  echo ""
fi

# ---
# Run GATK ASEReadCounter
# ---
echo "--- Running GATK ASEReadCounter ---"
echo "  Input BAM: ${INPUT_BAM}"
echo "  Variants VCF (het only): ${VARIANTS_HET_VCF}"
echo "  Reference FASTA: ${REFERENCE_FASTA}"
echo "  Sample: ${SAMPLE_NAME}"
echo "  Output TAB: ${FINAL_OUTPUT_NAME}"
echo "  Min Mapping Quality: ${MIN_MQ}"
echo "  Min Base Quality: ${MIN_BQ}"
echo "  Perform Indels: ${PERFORM_INDELS_FLAG:-No}"
echo "  Verbosity: ${VERBOSITY_LEVEL}"
echo "---"

gatk ASEReadCounter \
  -R "${REFERENCE_FASTA}" \
  -V "${VARIANTS_HET_VCF}" \
  -I "${INPUT_BAM}" \
  -sample "${SAMPLE_NAME}" \
  -O "${FINAL_OUTPUT_NAME}" \
  --min-mapping-quality "${MIN_MQ}" \
  --min-base-quality "${MIN_BQ}" \
  ${PERFORM_INDELS_FLAG} \
  --verbosity "${VERBOSITY_LEVEL}" || { echo "Error: ASEReadCounter failed."; exit 1; }

echo "---"
echo "ASEReadCounter completed successfully."
echo "Output saved to: ${FINAL_OUTPUT_NAME}"
echo "---"
