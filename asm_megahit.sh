#!/bin/bash
set -o errexit;

# Co-assembly of metagenome - MEGAHIT
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Description
# This script will run MEGAHIT using all available paired-end and single-end libraries by merging each together.

# Usage
# asm_megahit.sh <trim_fq_dir> <out_dir> <threads> 

# Dependencies
# - MEGAHIT
# - GNU parallel

# Citations

# Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W.
# MEGAHIT v1.0: A fast and scalable metagenome assembler driven by advanced methodologies and community practices.
# https://doi.org/10.1016/j.ymeth. (2016)

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
fq_dir=$1
out_dir=$2
threads=$3


# Merge quality-trimmed paired-end and single-end libraries
parallel "cat $fq_dir/*_[rR]{}* > $out_dir/merged_R{}.fq.gz" ::: 1 2
cat $fq_dir/*singles* > $out_dir/merged_singles.fq.gz

# De novo metagenome assembly
megahit --presets meta-sensitive -t $threads -1 $out_dir/merged_R1.fq.gz -2 $out_dir/merged_R2.fq.gz -r $out_dir/merged_singles.fq.gz -o $out_dir/asm/

# Organize files
mv $out_dir/asm/* $out_dir
mv $out_dir/final.contigs.fa $out_dir/asm_contigs_megahit.fasta

rm -r $out_dir/asm/
rm $out_dir/*fq.gz
