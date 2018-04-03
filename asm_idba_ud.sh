#!/bin/bash
set -o errexit;

# Co-assembly of metagenome - IDBA_UD
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Description
# This script will run IDBA_UD using all available paired-end libraries by merging each together.

# Usage
# asm_idba_ud.sh <trim_fq_dir> <out_dir> <threads> 

# Dependencies
# - IDBA_UD
# - GNU parallel

# Citations

# Yu Peng, Henry C. M. Leung, S. M. Yiu, Francis Y. L. Chin; 
# IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth. 
# Bioinformatics (2012)

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
fq_dir=$1
out_dir=$2
threads=$3


# Merge quality-trimmed paired-end libraries
parallel "cat $fq_dir/*_[rR]{}* > $out_dir/merged_R{}.fq.gz" ::: 1 2

# Convert to interleaved FASTA format
fq2fa --merge <(gunzip -c $out_dir/merged_R1.fq.gz) <(gunzip -c $out_dir/merged_R2.fq.gz) $out_dir/merged_pe.fa
rm $out_dir/*fq.gz

# De novo metagenome assembly
idba_ud --pre_correction --mink 33 --maxk 213 --step 20 --num_threads $threads --long_read $out_dir/merged_pe.fa --out $out_dir

# Organize files
mv $out_dir/contig.fa $out_dir/asm_contigs_idba_ud.fasta
rm $out_dir/merged_pe.fa
