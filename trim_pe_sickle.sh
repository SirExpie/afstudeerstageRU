#!/bin/bash
set -o errexit;

# Quality based trimming of paired-end reads - Sickle
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# trim_pe_sickle.sh <fq_dir> <q_treshold> <min_ctg_length> <out_dir>

# Dependencies
# - Sickle
# - GNU parallel
# - Perl 'rename' function

# Citations

# Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files 
# (Version 1.33) [Software]. Available at https://github.com/najoshi/sickle.

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
fq_dir=$1
q_tresh=$2
min_len=$3
out_dir=$4


# Perform trimming
parallel --link "sickle pe -t sanger -q $q_tresh -l $min_len -n -g --quiet -f {1} -r {2} -o $out_dir/{=1 s:.*/::;s:\..*$:_trim.fq.gz:; =} -p $out_dir/{=2 s:.*/::;s:\..*$:_trim.fq.gz:; =} -s $out_dir/{=1 s:.*/::;s:_[rR]1::;s:\..*$:_singles_trim.fq.gz:; =}" ::: $(ls $fq_dir/*_[rR]1*) ::: $(ls $fq_dir/*_[rR]2*)
