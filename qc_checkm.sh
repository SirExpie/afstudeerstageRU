#!/bin/bash
set -o errexit;

# Assessing quality of genome bins - CheckM
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# qc_checkm.sh <prefix> <contigs.fa> <bin_dir> <extension> <map_dir> <out_dir> <threads>
# - <bin_dir> directory that contains all bins: *.fasta files
# - <extension> extension of fasta file (.fa, .fna, .fasta, etc.)
# - <map_dir> directory that contains sorted, indexed bam files
# # Note: <out_dir> must be empty

# Dependencies
# - CheckM

# Citations

# Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P. & Tyson, G. W. 
# CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. 
# Genome research 25, 1043-1055, doi:10.1101/gr.186072.114 (2015).

# ------------------------------------------------------------------------------

prefix=$1
asm=$2
bin_dir=$3
extension=$4
map_dir=$5
out_dir=$6
threads=$7


# Activate virtual environment
source activate binmate

# Place bins in the reference genome tree
checkm tree -q -t $threads --pplacer_threads $threads -x $extension $bin_dir $out_dir

# Assess phylogenetic markers found in each bin
checkm tree_qa -o 2 --tab_table -q $out_dir > $out_dir/$prefix.phylo_markers.tsv

# Infer lineage-specific marker sets for each bin
checkm lineage_set -q $out_dir $out_dir/markers.mf

# Identify marker genes in bins
checkm analyze -q -t $threads $out_dir/markers.mf -x $extension $bin_dir $out_dir

# Calculate coverage of sequences
checkm coverage --all_reads --min_align 0.95 -q -t $threads -q -x $extension $bin_dir $out_dir/$prefix.coverage.tsv $map_dir/*.sorted.bam

# Calculate percentage of reads mapped to each bin
checkm profile --tab_table -q $out_dir/$prefix.coverage.tsv > $out_dir/$prefix.mapped_reads.tsv

# Assess bins for contamination and completeness
checkm qa -o 2 -t $threads -q --tab_table -c $out_dir/$prefix.coverage.tsv $out_dir/markers.mf $out_dir > $out_dir/$prefix.quality.tsv

# Join tab-separated tables
checkm join_tables $out_dir/$prefix.phylo_markers.tsv $out_dir/$prefix.quality.tsv $out_dir/$prefix.mapped_reads.tsv > $out_dir/checkm_results.$prefix.tsv

# Bar plot of bin completeness, contamination, and strain heterogeneity
checkm bin_qa_plot --image_type pdf --dpi 300 -q  -x $extension $out_dir $bin_dir $out_dir
mv $out_dir/bin_qa_plot.pdf $out_dir/$prefix.bin_qa_plot.pdf

# Identify unbinned contigs
checkm unbinned -q -x $extension $bin_dir $asm $out_dir/$prefix.unbinned_contigs.fasta $out_dir/$prefix.unbinned_stats.tsv

# Identify SSU (16S/18S) rRNAs in sequences
checkm ssu_finder -q -t $threads -x $extension $asm $bin_dir $out_dir/$prefix_16s_seq/

# Deactivate virtual environment
source deactivate
