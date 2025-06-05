#!/bin/bash

# ---
# Script to add or replace Read Group (RG) information in a BAM file
# using GATK's AddOrReplaceReadGroups.
#
# Requirements:
# - GATK (Genome Analysis Toolkit) must be in your system's PATH.
# ---

# Function to display usage information
usage() {
  echo "Usage: $0 <input_bam> <rgsm_value> [options]"
  echo ""
  echo "Arguments:"
  echo "  <input_bam>        Path to the input BAM file (e.g., /path/to/reads.bam)"
  echo "  <rgsm_value>       REQUIRED: The biological sample name (SM tag)."
  echo "                     This MUST match the sample column name in your VCF file."
  echo ""
  echo "Options (for Read Group tags):"
  echo "  --rgid <ID>        Read Group ID (ID tag). Default: base name of input BAM."
  echo "  --rglb <LIB>       Read Group Library (LB tag). Default: base name of input BAM."
  echo "  --rgpl <PLATFORM>  Read Group Platform (PL tag). Default: 'Illumina'."
  echo "  --rgpu <UNIT>      Read Group Platform Unit (PU tag). Default: base name of input BAM + '.unit'."
  echo "  -o, --output <file> Specify output BAM file name. Default: <input_bam_base>.RG.bam."
  echo "  -h, --help         Display this help message."
  echo ""
  echo "Example:"
  echo "  $0 my_raw_data.bam hESC_H9"
  echo "  $0 reads.bam SampleX --rglb lib_prep1 --rgpl 'DNBSEQ' --rgpu 'flowcellA.lane1'"
  exit 1
}

# ---
# Parse command-line arguments
# ---

# Check for help flag first
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
fi

# Assign required positional arguments
INPUT_BAM="$1"
RGSM="$2" # REQUIRED and described first as per request

# Shift arguments to process options
shift 2

# Validate required positional arguments
if [[ -z "$INPUT_BAM" || -z "$RGSM" ]]; then
  echo " "
  echo "Error: Missing required arguments. Please provide input BAM and RGSM value."
  usage
fi

# Validate input BAM file
if [[ ! -f "$INPUT_BAM" ]]; then
  echo " "
  echo "Error: Input BAM file '${INPUT_BAM}' does not exist or is not a regular file."
  exit 1
fi

# Derive base filename and directory for defaults
BASE_FILENAME=$(basename "${INPUT_BAM%.bam}") # e.g., "my_raw_data" from "/path/to/my_raw_data.bam"
INPUT_BAM_DIR=$(dirname "$INPUT_BAM")         # e.g., "/path/to"

# Set default values for Read Group tags and output file
RGID="${BASE_FILENAME}"                               # Default: base name of input BAM
RGLB="${BASE_FILENAME}"                               # Default: base name of input BAM
RGPL="Illumina"                                       # Default: Illumina
RGPU="${BASE_FILENAME}.unit"              	      # Default: base name + .unit
OUTPUT_BAM="${INPUT_BAM_DIR}/${BASE_FILENAME}.RG.bam" # Default output file

# Parse optional arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --rgid)
      RGID="$2"
      shift 2
      ;;
    --rglb)
      RGLB="$2"
      shift 2
      ;;
    --rgpl)
      RGPL="$2"
      shift 2
      ;;
    --rgpu)
      RGPU="$2"
      shift 2
      ;;
    -o|--output)
      OUTPUT_BAM="$2"
      shift 2
      ;;
    *)
      echo "Error: Unknown option '$1'"
      usage
      ;;
  esac
done

# ---
# Run GATK AddOrReplaceReadGroups
# ---

echo "--- Read Group Details ---"
echo "  Input BAM:                     ${INPUT_BAM}"
echo "  Output BAM:                    ${OUTPUT_BAM}"
echo "  Read Group ID (ID):            ${RGID}"
echo "  Read Group Library (LB):       ${RGLB}"
echo "  Read Group Platform (PL):      ${RGPL}"
echo "  Read Group Platform Unit (PU): ${RGPU}"
echo "  Read Group Sample Name (SM):   ${RGSM}" # This remains the primary focus
echo "---"

# Construct the GATK command
GATK_COMMAND="gatk AddOrReplaceReadGroups"
GATK_COMMAND+=" -I \"${INPUT_BAM}\""
GATK_COMMAND+=" -O \"${OUTPUT_BAM}\""
GATK_COMMAND+=" --RGID \"${RGID}\""
GATK_COMMAND+=" --RGLB \"${RGLB}\""
GATK_COMMAND+=" --RGPL \"${RGPL}\""
GATK_COMMAND+=" --RGPU \"${RGPU}\""
GATK_COMMAND+=" --RGSM \"${RGSM}\""
# Removed RGCN, RGDT, RGDS, RGPG parameters
GATK_COMMAND+=" --SORT_ORDER coordinate"
GATK_COMMAND+=" --CREATE_INDEX true"

# Execute the GATK command
echo "Executing: ${GATK_COMMAND}"
eval "${GATK_COMMAND}" || \
  { echo "Error: GATK AddOrReplaceReadGroups failed. Please check GATK output above."; exit 1; }

echo "---"
echo "Read groups added successfully."
echo "Output saved to: ${OUTPUT_BAM}"
echo "Index created: ${OUTPUT_BAM}.bai"
echo "---"
