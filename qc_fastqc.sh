#!/bin/bash
set -o errexit;

# Sequencing data quality control - FastQC
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# qc_fastqc.sh <fq_dir> <out_dir>

# Dependencies
# - FastQC
# - GNU parallel
# - Java Runtime Environment (JRE, readily available on most systems)

# Citations

# Andrews S. (2010). 
# FastQC: a quality control tool for high throughput sequence data. 
# Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
fq_dir=$1
out_dir=$2


parallel "fastqc -f fastq --extract --quiet -o $out_dir {}" ::: $(find -L $fq_dir -maxdepth 1 -type f -regextype posix-extended -iregex '.*/.*\.fastq.*|.*/.*\.fq.*')
