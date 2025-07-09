# Allele Specific Expression (ASE) from Stranded RNA-seq

Runs a simple pipeline for performing Allele-Specific Expression (ASE) analysis on data derived from stranded RNA-seq. Takes a BAM file, a VCF file, and an annotation GTF file, and produces an output table containing SNP-level allele counts (reference and alternate), genomic coordinates for each heterozygous SNP, allelic imbalance metrics, and computed minor allele frequencies (MAF). This output enables downstream analyses such as detection of allele-specific expression patterns and gene-level aggregation of ASE signals.

## Install Required Packages (Conda Env)

To create a new conda environment named "GATK4" containing all required packages run: 

`bash /src/run-installPackages.sh`

This Bash script creates a Conda environment named gatk4, configures the necessary Conda channels with strict priority, activates the environment, installs GATK4, and then installs additional tools (bedtools, samtools, gawk, coreutils, pandas, numpy, matplotlib) required for running various bioinformatics scripts.

## Quick Start: Run the stranded-ASE-analysis pipeline

This is the main pipeline script for ASE analysis from a single stranded RNA-seq BAM file.

`bash run-strandedASE.sh <input.RG.bam> <variants.vcf> <genome.fa> <annotation.gtf> [options]`

This script automates all bioinformatics steps (src folder) into a single command:

1. splitBAM: Divides the input BAM into forward and reverse strand-specific BAM files.
2. aseReadCounter: Runs ASEReadCounter on both strands against variants on the VCF file.
3. annotateVariants: Suplements ASEReadCounter ouput with gene information
4. geneLevelMAF: Aggregates annotated SNP-level data to calculate gene-level MAF.

**Required Arguments**: <br>
- <input_bam> : Path to the input BAM file (must contain Read Groups!). <br>
- <variants_vcf> : Path to the VCF file with variant sites. <br>
- <reference_fasta> : Path to the reference FASTA file. <br>
- <reference_gtf> : Path to the reference GTF annotation file. <br>

Optional Arguments: <br>
These are passed directly to the run-aseReadCounter.sh sub-script: <br>
- --min-mapping-quality <INT>: Minimum mapping quality for reads (default: 20).  <br>
- --min-base-quality <INT>: Minimum base quality for bases (default: 20).  <br>
- --perform-indels: Flag to include indels in allele counting.  <br>
- --verbosity <LEVEL>: GATK verbosity level (default: INFO).  <br>
- -h, --help: Display usage information.  <br>

**MAF-Based Metrics in the Output Table**:

- **minor_fq**: Mean of SNP-level MAFs per gene

For each SNP in a gene, computes `MAF_SNP = min(refCount, altCount) / totalCount`. <br>
These values are then averaged across all SNPs in that gene. This represents the average allelic balance across sites.

- **minor_fq_tot**: Gene-level MAF based on aggregated counts

First, all refCount and altCount are summed across SNPs for each gene. <br> 
Then computes `minor_fq_tot = min(total_refCount, total_altCount) / totalCount_sum`. <br>
This captures the overall balance across the gene (possibly masking outlier SNPs).

**Expression Classification:**

The expression status of each gene is assigned using `minor_fq`: <br>
- biallelic if MAF ≥ 0.40 <br>
- imbalance if 0.10 < MAF < 0.40 <br>
- monoallelic if MAF ≤ 0.10 <br>


## Pre-Processing Step: Add Read Groups to BAM Files (if missing)
GATK4 tools require properly formatted read groups in your BAM files to function correctly and ensure accurate downstream analysis (e.g., duplicate marking, variant calling, quality control).
If your bam file does not contain read groups, you will get an error message. To verify if your BAM file already contains read group information, use the following command:

`samtools view -H your_input.bam | grep '^@RG'` <br>

If this command produces no output, or if the output is incomplete or incorrect, you need to add or replace read groups. <br>
If it shows well-defined @RG lines, you can skip this step. <br>

*Run the wrapper for GATK's AddOrReplaceReadGroups.*

`bash /src/run-addReadGroup.sh path/to/your_input.bam --rgsm_value 'sample_name'`

The only required argument here is the rgsm_value (sample name) that will be added to the Read Group.  <br>
It is critical because it MUST match the sample column name in your VCF file for subsequent analysis steps. <br>
This script will create a new BAM file (filename.RG.bam) in the same directory as your input file. <br>

Optional Arguments: <br>
--rgid <ID>: Read Group ID (ID tag). Default: base name of input BAM. <br>
--rglb <LIB>: Read Group Library (LB tag). Default: base name of input BAM. <br> 
--rgpl <PLATFORM>: Read Group Platform (PL tag). Default: Illumina. Valid values include ILLUMINA, SOLID, PACBIO, etc. <br>
--rgpu <UNIT>: Read Group Platform Unit (PU tag). Default: base name of input BAM + .unit. <br>
 -o <file>: Specify the output BAM file name. Default: <input_bam_base>.RG.bam in the input BAM's directory. <br>


## Background Functions

### 1. splitBAM
Divides the input BAM into forward and reverse strand-specific BAM files.

GATK's ASEReadCounter doesn't inherently account for strand-specific read orientation. Our workaround involves splitting the aligned BAM file into forward and reverse strand components, then running ASEReadCounter individually on each. To run this step, simply type:

`bash /src/run-splitBAM.sh <input_bam>`

Output:
- Creates two new bam files (fwd and rev) inside a new folder (default: /data/bams-split)

Optional Arguments: <br>
- --threads <INT>: Number of cores to be used (default: 8).  <br>
- --output-dir <INT>: Specify output folder name (default: /data/bams-split).  <br>


### 2. aseReadCounter
Runs ASEReadCounter on both strands against variants on the VCF file.

This step calculates allele-specific read counts for the variants of interest from your prepared BAM files. Use the run-aseReadCounter.sh script with the following arguments: 

`bash /src/run-aseReadCounter.sh <input_bam> <variants.vcf.gz> <reference.fasta>`

Output:
- A new folder named out-aseReadCounter is created in the directory where the script is executed.
- Inside, a *tab-delimited file* is created containing **allele-specific read counts for each variant position**.

### 3. annotateVariants
Suplements ASEReadCounter ouput with gene information

This step annotates Allele-Specific Expression (ASE) count data—produced by GATK’s ASEReadCounter—with gene information. It takes as input two allele count files (in .tab format) corresponding to the forward and reverse strands, along with a gene annotation file in .gtf format.

`bash /src/run-annotateVariants.sh <fwd_tab_file> <rev_tab_file> <gtf_file>`

Output:
- Produces a .tab file named `filename.annotated.tab` combining ASE counts with gene annotations. 

**Key Steps in the Script:**

1. Converts the gtf annotation and the ASEReadCounter output files into BED format.
2. Splits the annotation BED into forward (fwd) and reverse (rev) strand-specific files for strand-aware analysis.
3. Uses *bedtools intersect -wa -wb* to map allelic counts to gene annotations in a strand-specific manner.
4. Creates .annotated.tab files (fwd and rev) combining ASE counts with gene annotations.
5. Both files are merged in the same table and saved in the same directory as the input files.

### 4. geneLevelMAF
Computes gene-level Minor Allele Frequency (MAF) and use it as a proxy to infer allelic imbalance. 

`python src/run-geneLevelMAF.py filename.annotated.tab`

**Optional arguments**: <br>
- --threshold: Minimum total counts for a SNP to be included. SNPs with 'totalCount' <= threshold will be excluded. Default: 10.
- --output: Path to output file. Default: '<input_filename>.geneLevelMAF'.

**Key Steps in the Script**:
1. Read SNP-level ASE data.
2. Filter SNPs with low total read counts.
3. Calculate per-SNP MAF and compute `minor_fq`.
4. Sum SNP counts per gene and compute `minor_fq_tot`
5. Concatenate all SNPs variant considered into a column.
6. Classify genes by expression pattern besed on the `minor_fq`.

## Helper Functions

#### Filter VCF file for heterozygous variants

`bash src/run-filterVCF.sh <filename.vcf.gz> <sample_name>`

- It uses GATK SelectVariants and writes the output to filename.HET.vcf.gz
- ASEReadCounter runs faster if only heterozygous positions are provided

#### Correcting Mapping Bias with WASP

When aligning DNA reads to a reference genome, differences between the sample and reference can cause reference allele bias, where reads matching the reference are more likely to align. This can distort allele-specific expression analyses. STAR includes an efficient implementation of WASP to correct for this bias when sample genotypes are provided (`--varVCFfile`). To use it, enable `--waspOutputMode` along with `--varVCFfile`.

```
STAR \
--runThreadN 8 \
--genomeDir $index \
--readFilesCommand zcat \
--readFilesIn file_1.fq.gz file_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1 \
--quantMode GeneCounts \
--sjdbGTFfile $annot_file \
--waspOutputMode SAMtag \
--varVCFfile $vcf_file \
--outSAMattributes All
```

Next, we can use `samtools` to select reads that are uniquely mapped (`-q 255`) and have allele-specific bias corrected — either they don’t overlap SNPs (`![vW]`) or passed the WASP filter (`[vW]==1`):

`samtools view -h -b -e '![vW] || [vW]==1' -q 255 -o filename.WASPfilter.aligned.bam filename.aligned.bam`

Overall, it’s safe to enable `--waspOutputMode` by default, then the vW tag to filter reads only when reference bias matters (e.g., in ASE or imprinting studies); otherwise, you can ignore it and use all reads normally.
But this procedure increases runtime and is only necessary for allele-specific analyses. 
