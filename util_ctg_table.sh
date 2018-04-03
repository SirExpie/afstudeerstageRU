#!/bin/bash
set -o errexit;

# Generate contig information table
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# util_ctg_table.sh <contigs.fasta> <map_dir> <bin_dir> <out_dir> <threads>

# Dependencies
# - Kaiju (with suitable index)
# - MetaBAT (jgi_summarize_bam_contig_depths script)

# Citations

# Menzel, P., Ng, K. L. & Krogh, A. 
# Fast and sensitive taxonomic classification for metagenomics with Kaiju. 
# Nat Commun 7, 11257, doi:10.1038/ncomms11257 (2016).

# Kang, D. D., Froula, J., Egan, R. & Wang, Z.
# MetaBAT, an efficient tool for accurately reconstructing single genomes from complex microbial communities.
# PeerJ 3, e1165, doi:10.7717/peerj.1165 (2015).

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
asm=$1
map_dir=$2
bin_dir=$3
out_dir=$4
threads=$5


prefix=$(basename $asm | sed 's/\.fa.*//')
kaiju_db=/data/data/databases/kaiju

# Generate contig data
$(dirname $(readlink -f $0))/qc_asm_stats.py -i $asm --datatable -o $out_dir/dat_$prefix.tsv

# Get contig bin assignments
$(dirname $(readlink -f $0))/util_bin2tsv.sh $bin_dir $out_dir/bin_$prefix.tsv

# Generate coverage data
jgi_summarize_bam_contig_depths --outputDepth $out_dir/cov_$prefix.tsv $map_dir/*.sorted.bam

# Classify contigs
kaiju -a greedy -e 5 -x -z $threads -t $kaiju_db/nodes.dmp -f $kaiju_db/kaiju_db_nr_euk.fmi -i $asm -o $out_dir/kaiju.out
addTaxonNames -p -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp -i $out_dir/kaiju.out -o $out_dir/tax_$prefix-taxonpath.tsv
addTaxonNames -r superkingdom,phylum,order,class,family,genus,species -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp -i $out_dir/kaiju.out -o $out_dir/tax_$prefix-split.tsv

sed -i 's/; /\t/g' $out_dir/tax_$prefix-split.tsv
rm $out_dir/kaiju.out

# R script

# rm $out_dir/dat_$prefix.tsv $out_dir/bin_$prefix.tsv $out_dir/cov_$prefix.tsv $out_dir/tax_$prefix-taxonpath.tsv $out_dir/tax_$prefix-spit.tsv


