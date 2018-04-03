#!/bin/bash
set -o errexit;

# Metagenome binning - MetaBAT2
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# bin_metabat2.sh <contigs.fasta> <min_ctg_len> <map_dir> <out_dir> <threads>

# Dependencies
# - MetaBAT2
# - Anaconda with "binmate" virtual environment

# Citations

# Kang, D. D., Froula, J., Egan, R. & Wang, Z. 
# MetaBAT, an efficient tool for accurately reconstructing single genomes from complex microbial communities. 
# PeerJ 3, e1165, doi:10.7717/peerj.1165 (2015).

# Anaconda Software Distribution. Computer software. Continuum Analytics. Web. <https://continuum.io>.

# ------------------------------------------------------------------------------

# Arguments
asm=$1
min_ctg_len=$2
map_dir=$3
out_dir=$4
threads=$5


mkdir -p $out_dir/bins_metabat2

# Activate virtual environment
source activate binmate

# Generage coverage profiles
jgi_summarize_bam_contig_depths --outputDepth $out_dir/depth_metabat2.tsv --pairedContigs $out_dir/paired_metabat2.tsv --minContigDepth 2 $map_dir/*.sorted.bam

# Run MetaBAT
metabat2 -a $out_dir/depth_metabat2.tsv -t $threads -i $asm -m $min_ctg_len -o $out_dir/bins_metabat2/metabat2

# Deactivate virtual environment
source deactivate
