#!/bin/bash

# ---
# Filters a VCF for PASS heterozygous SNPs only for a given sample name.
# Uses GATK SelectVariants and writes the output to <input>.HET.SNPs.vcf.gz.
#
# Usage: ./run-filterVCF.sh <VCF.gz> <SAMPLE_NAME>
# Example: ./run-filterVCF.sh cohort.vcf.gz Sample_01
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
    echo "This script filters for PASS heterozygous SNPs only for the specified sample"
    echo "using GATK SelectVariants. It creates an output file named"
    echo "<VCF>.HET.SNPs.vcf.gz in the same directory."
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
VARIANTS_HET_VCF="${VARIANTS_VCF/.vcf.gz/.HET.SNPs.vcf.gz}"

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

# Filter for PASS heterozygous SNPs only
if [ ! -f "$VARIANTS_HET_VCF" ]; then
    echo "🔧 Creating PASS heterozygous SNP-only VCF for sample '$SAMPLE_NAME': $VARIANTS_HET_VCF"

    gatk SelectVariants \
        -V "$VARIANTS_VCF" \
        --select-type-to-include SNP \
        --exclude-filtered \
        --select "vc.getGenotype('${SAMPLE_NAME}').isHet()" \
        -O "$VARIANTS_HET_VCF"

    echo "✅ Done: $VARIANTS_HET_VCF"
else
    echo "⚠️  File already exists: $VARIANTS_HET_VCF. Skipping generation."
fi