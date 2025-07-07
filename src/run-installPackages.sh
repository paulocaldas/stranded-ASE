#!/bin/bash
# This script creates a new conda environment and installs all necessary packages

# Configure Conda Channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# This line is crucial for 'conda activate' to work in scripts
source "$(conda info --base)/etc/profile.d/conda.sh"

# Create a conda environment (initially empty)
conda create -n gatk4 -y

# Activate the environment
conda activate gatk4

# Install GATK4 in a separate step
conda install -y gatk4

# Install all other required packages for all scripts
conda install -y bedtools bcftools samtools gawk coreutils pandas numpy matplotlib
