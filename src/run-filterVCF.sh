#!/bin/bash

# ---
# This script filters a VCF file for heterozygous variants for a given sample name.
# It uses GATK SelectVariants and writes output to <input>.HET.vcf.gz.
# Usage: ./filter_het_variants.sh <VCF.gz> <SAMPLE_NAME>
# ---

set -euo pipefail

VARIANTS_VCF="$1"
SAMPLE_NAME="$2"
VARIANTS_HET_VCF="${VARIANTS_VCF/.vcf.gz/.HET.vcf.gz}"

# Check input
if [ ! -f "$VARIANTS_VCF" ]; then
    echo "❌ ERROR: VCF file not found: $VARIANTS_VCF"
    exit 1
fi

# Check if sample exists in VCF
if ! bcftools query -l "$VARIANTS_VCF" | grep -q "^${SAMPLE_NAME}$"; then
    echo "❌ ERROR: Sample '${SAMPLE_NAME}' not found in VCF."
    exit 1
fi

# Filter heterozygous variants
if [ ! -f "$VARIANTS_HET_VCF" ]; then
    echo "🔧 Creating heterozygous-only VCF: $VARIANTS_HET_VCF"

    gatk SelectVariants \
        -V "$VARIANTS_VCF" \
        --select "vc.getGenotype('${SAMPLE_NAME}').isHet()" \
        -O "$VARIANTS_HET_VCF"

    echo "✅ Done: $VARIANTS_HET_VCF"
else
    echo "⚠️  File already exists: $VARIANTS_HET_VCF, skipping."
fi
