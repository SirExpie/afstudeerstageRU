#!/bin/bash
set -o errexit;

# Metagenome binning - MetaBAT
# ------------------------------------------------------------------------------
# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# bin_binsanity.sh <contigs.fasta> <map_dir> <threads>

# Dependencies
# - BinSanity
# - CheckM
# - Anaconda with "binmate" virtual environment

# Note
# BinSanity writes results to the present working directory (even when using -o option)

# Citations

# Graham, E. D., Heidelberg, J. F. & Tully, B. J. 
# BinSanity: unsupervised clustering of environmental microbial assemblies using coverage and affinity propagation. 
# PeerJ 5, e3035, doi:10.7717/peerj.3035 (2017).

# Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P. & Tyson, G. W. 
# CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. 
# Genome research 25, 1043-1055, doi:10.1101/gr.186072.114 (2015).

# Anaconda Software Distribution. Computer software. Continuum Analytics. Web. <https://continuum.io>.

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
asm=$1
map_dir=$2 #map/pe/
threads=$3

# Activate virtual environment
source activate binmate

# Extract contig IDs
echo $asm
echo "dirname: "$(dirname $asm)" basename: "$(basename $asm)
get-ids -f $(dirname $asm) -l $(basename $asm) -o contig_ids.txt

# Generate coverage profiles
# Note: if extension is not *.bam, Binsanity-profile/Subread will fail
# The "scale" transformation of coverage values is recommended
Binsanity-profile -i $asm --ids contig_ids.txt -s $map_dir --transform scale -c cov_binsanity

# Run BinSanity
Binsanity-wf -f $(dirname $asm) -l $(basename $asm) -c *.x100.lognorm --threads $threads

# Move and rename bins
echo "Moving, renaming folders..."
mv BINSANITY-RESULTS/BinSanity-Final-bins/ .
mv BinSanity-Final-bins/ bins_binsanity/

cd bins_binsanity/
parallel 'mv {} binsanity.{.}.fa' ::: $(ls)

# Deactivate virtual environment
source deactivate
