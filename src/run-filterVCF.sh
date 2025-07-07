#!/bin/bash

# ---
# This script filters a VCF file for heterozygous variants for a given sample name.
# It uses GATK SelectVariants and writes the output to <input>.HET.vcf.gz.
#
# Usage: ./filter_het_variants.sh <VCF.gz> <SAMPLE_NAME>
# Example: ./filter_het_variants.sh cohort.vcf.gz Sample_01
# ---

set -euo pipefail

# Function to show usage help
usage() {
    echo ""
    echo "Usage: $0 <VCF.gz> <SAMPLE_NAME>"
    echo ""
    echo "Arguments:"
    echo "  <VCF.gz>       Path to the input compressed VCF file (.vcf.gz)"
    echo "  <SAMPLE_NAME>  Sample name present in the VCF header (FORMAT columns)"
    echo ""
    echo "This script filters heterozygous variants for the specified sample using GATK SelectVariants."
    echo "It creates an output file named <VCF>.HET.vcf.gz in the same directory."
    echo ""
    exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "❌ ERROR: Incorrect number of arguments."
    usage
fi

VARIANTS_VCF="$1"
SAMPLE_NAME="$2"
VARIANTS_HET_VCF="${VARIANTS_VCF/.vcf.gz/.HET.vcf.gz}"

# Check if input VCF exists
if [ ! -f "$VARIANTS_VCF" ]; then
    echo "❌ ERROR: VCF file not found: $VARIANTS_VCF"
    exit 1
fi

# Check if input VCF is indexed
if [ ! -f "${VARIANTS_VCF}.tbi" ]; then
    echo "❌ ERROR: Index (.tbi) file not found for VCF. Please index with 'bcftools index $VARIANTS_VCF'"
    exit 1
fi

# Check if sample exists in VCF
if ! bcftools query -l "$VARIANTS_VCF" | grep -qx "${SAMPLE_NAME}"; then
    echo "❌ ERROR: Sample '${SAMPLE_NAME}' not found in VCF."
    echo "👉 Available samples in VCF:"
    bcftools query -l "$VARIANTS_VCF"
    exit 1
fi

# Filter heterozygous variants for the given sample
if [ ! -f "$VARIANTS_HET_VCF" ]; then
    echo "🔧 Creating heterozygous-only VCF for sample '$SAMPLE_NAME': $VARIANTS_HET_VCF"

    gatk SelectVariants \
        -V "$VARIANTS_VCF" \
        --select "vc.getGenotype('${SAMPLE_NAME}').isHet()" \
        -O "$VARIANTS_HET_VCF"

    echo "✅ Done: $VARIANTS_HET_VCF"
else
    echo "⚠️  File already exists: $VARIANTS_HET_VCF. Skipping generation."
fi
