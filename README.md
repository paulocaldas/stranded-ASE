# Allele Specific Expression (ASE) from Stranded RNA-seq

Runs a simple pipeline for performing Allele Specific Expression (ASE) analysis on data derived from stranded RNA-seq

## Install Required Packages (Conda Env)

To create a new conda environment named "GATK4" containing all required packages run: 

`./src/run-installPackages.sh`

## Quick Start: Run the stranded-ASE-analysis pipeline

This is the main pipeline script for ASE analysis from a single stranded RNA-seq BAM file.

`./run-strandedASE.sh file path/to/your_input.RG.bam variants.vcf.gz genome.fa genome.gtf`

This script automates all bioinformatics steps (src folder) into a single command:

1. splitBAM: Divides the input BAM into forward and reverse strand-specific BAM files.
2. aseReadCounter: Runs ASEReadCounter on both strands against variants on the VCF file.
3. annotateVariants: Suplements ASEReadCounter ouput with gene information
4. geneLevelMAF: Aggregates annotated SNP-level data to calculate gene-level MAF.

Required Arguments: <br>
**<input_bam>**: Path to the input BAM file (must contain Read Groups). <br>
**<variants_vcf>**: Path to the VCF file with variant sites. <br>
**<reference_fasta>**: Path to the reference FASTA file. <br>
**<reference_gtf>**: Path to the reference GTF annotation file. <br>

Optional Arguments: <br>
These are passed directly to the run-aseReadCounter.sh sub-script: <br>
--min-mapping-quality <INT>: Minimum mapping quality for reads (default: 20).  <br>
--min-base-quality <INT>: Minimum base quality for bases (default: 20).  <br>
--perform-indels: Flag to include indels in allele counting.  <br>
--verbosity <LEVEL>: GATK verbosity level (default: INFO).  <br>
-h, --help: Display usage information.  <br>

## Pre-Processing Step: Add Read Groups to BAM Files (if missing)
GATK4 tools require properly formatted read groups in your BAM files to function correctly and ensure accurate downstream analysis (e.g., duplicate marking, variant calling, quality control).
If your bam file does not contain read groups, you will get an error message. To verify if your BAM file already contains read group information, use the following command:

`samtools view -H your_input.bam | grep '^@RG'` <br>

If this command produces no output, or if the output is incomplete or incorrect, you need to add or replace read groups. <br>
If it shows well-defined @RG lines, you can skip this step. <br>

*Run the wrapper for GATK's AddOrReplaceReadGroups.*

`./src/run-addReadGroup.sh path/to/your_input.bam --rgsm_value 'sample_name'`

The only required argument here is the rgsm_value (sample name) that will be added to the Read Group.  <br>
It is critical because it MUST match the sample column name in your VCF file for subsequent analysis steps. <br>
This script will create a new BAM file (filename.RG.bam) in the same directory as your input file. <br>

Optional Arguments: <br>
--rgid <ID>: Read Group ID (ID tag). Default: base name of input BAM. <br>
--rglb <LIB>: Read Group Library (LB tag). Default: base name of input BAM. <br> 
--rgpl <PLATFORM>: Read Group Platform (PL tag). Default: Illumina. Valid values include ILLUMINA, SOLID, LS454, HELICOS, PACBIO, etc. <br>
--rgpu <UNIT>: Read Group Platform Unit (PU tag). Default: base name of input BAM + .unit. <br>
 -o <file>: Specify the output BAM file name. Default: <input_bam_base>.RG.bam in the input BAM's directory. <br>



## Background Functions


