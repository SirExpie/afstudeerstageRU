#!/bin/bash
set -o errexit;

# Metagenome binning - CONCOCT
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# bin_concoct.sh <contigs.fasta> <map_dir> <out_dir> <threads>

# Dependencies
# - Anaconda with "binmate" virtual environment
# - CONCOCT
# - perl "rename" function

# Notes
# Make sure CONCOCT is installed within the virtual environment. 
# Otherwise will trigger "ImportError: No module named vbgmm".

# Citations

# Anaconda Software Distribution. Computer software. Continuum Analytics. Web. <https://continuum.io>.

# Alneberg, J. et al. 
# Binning metagenomic contigs by coverage and composition. 
# Nature methods 11, 1144-1146, doi:10.1038/nmeth.3103 (2014).

# ------------------------------------------------------------------------------

# Arguments
asm=$1
map_dir=$2
out_dir=$3
threads=$4

# Activate virtual environment
source activate binmate

mkdir -p $out_dir/bins_concoct

# The CONCOCT documentation recommends cutting up assembled contigs into 10 kb chunks
# This recommendation is ignored here

# Generate coverage profiles
bam_files=$(find $map_dir -name "*.sorted.bam")
gen_input_table.py $asm $bam_files | cut -f1,3- > $out_dir/cov_concoct.tsv

# Run CONCOCT
concoct --coverage_file $out_dir/cov_concoct.tsv --composition_file $asm -r 250 -b $out_dir/output/

# Extract clusters/bins
extract_fasta_bins.py --output_path $out_dir/bins_concoct/ $asm $out_dir/output/clustering_gt1000.csv

# Rename bins
cd $out_dir/bins_concoct/
rename 's/^/concoct\./' *.fa

# Generate and incorporate linkage information - ERROR: extraction of bins fails
# find $map_dir -name "*.sorted.bam" -exec basename {} \; > $out_dir/samples.txt
# bam_to_linkage.py -m $threads --regionlength 500 --fullsearch --samplenames $out_dir/samples.txt $asm $bam_files > $out_dir/linkage_concoct.tsv
# ClusterLinkNOverlap.pl --cfile=output/clustering_gt1000.csv --lfile=$out_dir/linkage_concoct.tsv --covfile=$out_dir/cov_concoct.tsv --quiet --ofile=output/clustering_gt1000_withLinkage.csv
# extract_fasta_bins.py --output_path bins ../../asm/contigs_gt2500.fasta concoct-output/clustering_gt1000_withLinkage.csv

# Deactivate virtual environment
source deactivate 
