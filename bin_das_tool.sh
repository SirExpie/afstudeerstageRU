#!/bin/bash
set -o errexit;

# Bin refinement - DAS Tool
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# bin_das_tool.sh <contigs.fasta> <data/bin_i.tsv,data/bin_n.tsv> <name_i,name_n> <out_dir> <threads>

# Dependencies
# - Anaconda with "binmate" virtual environment
# - DAS tool
# - GNU parallel

# Citations

# Sieber, M. K. et al. 
# Recovery of genomes from metagenomes via a dereplication, aggregation, and scoring strategy.
# bioRxiv preprint first posted online Feb. 11, 2017; doi: http://dx.doi.org/10.1101/107789.

# Anaconda Software Distribution. Computer software. Continuum Analytics. Web. <https://continuum.io>.

# O. Tange (2011): GNU Parallel - The Command-Line Power Too

# ------------------------------------------------------------------------------

# Arguments
asm=$1
bin_files=$2
bin_labels=$3
out_dir=$4
threads=$5


# Activate virtual environment
source activate binmate

# Run DAS Tool
DAS_Tool.sh --contigs $asm --bins $bin_files -l $bin_labels --write_bins 1 --threads $threads -o $out_dir/das_tool

# Rename bin folder + bins
cd bin/das_tool/
mv das_tool_DASTool_bins/ bins_das_tool/

cd bins_das_tool/
rename 's/\.contigs//' *fa
rename 's/(.*)$/das_tool\.$1/' *fa

# Deactivate virtual environment
source deactivate
