#!/bin/bash
set -o errexit;

# Taxonomic classification of sequences - Kaiju
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# classify_se_reads_kaiju.sh <prefix> <kaiju_db> <trim.fq.gz> <out_dir> <threads>

# Dependencies
# - Kaiju (with suitable index)
# - GNU parallel
# - Krona

# Citations

# Menzel, P., Ng, K. L. & Krogh, A. 
# Fast and sensitive taxonomic classification for metagenomics with Kaiju. 
# Nat Commun 7, 11257, doi:10.1038/ncomms11257 (2016).

# Ondov, B. D., Bergman, N. H. & Phillippy, A. M. 
# Interactive metagenomic visualization in a Web browser. 
# BMC bioinformatics 12, 385, doi:10.1186/1471-2105-12-385 (2011).

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
prefix=$1
kaiju_db=$2
r1=$3
out_dir=$4
threads=$5


# Classify reads
if [[ $r1 =~ \.gz$ ]]; then 
  kaiju -a greedy -e 5 -x -v -z $threads -t $kaiju_db/nodes.dmp -f $kaiju_db/kaiju_db_nr_euk.fmi -i <(gunzip -c $r1) -o $out_dir/$prefix.kaiju.out
else
  kaiju -a greedy -e 5 -x -v -z $threads -t $kaiju_db/nodes.dmp -f $kaiju_db/kaiju_db_nr_euk.fmi -i $r1 -o $out_dir/$prefix.kaiju.out
fi

addTaxonNames -p -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp -i $out_dir/$prefix.kaiju.out -o $out_dir/$prefix.classified-taxonpath.kaiju.txt
rm $out_dir/$prefix.kaiju.out

# Generate Krona interactive chart
kaiju2krona -u -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp -i $out_dir/$prefix.classified-taxonpath.kaiju.txt -o $out_dir/$prefix.kaiju.krona.out
ktImportText -o $out_dir/$prefix.krona.html $out_dir/$prefix.kaiju.krona.out
rm $out_dir/$prefix.kaiju.krona.out

# Generate summary reports
parallel "kaijuReport -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp -i $out_dir/$prefix.classified-taxonpath.kaiju.txt -r {} -o $out_dir/$prefix.summary-{}.kaiju.txt" ::: phylum class order family genus species
