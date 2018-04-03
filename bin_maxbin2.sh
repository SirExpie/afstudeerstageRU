#!/bin/bash
set -o errexit;

# Metagenome binning - MaxBin 2.0
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# bin_maxbin2.sh <contigs.fasta> <map_dir> <out_dir> <threads>

# Dependencies
# - Anaconda with "binmate" virtual environment
# - MaxBin 2.0
# - MetaBat (script: jgi_summarize_bam_contig_depths)
# - GNU parallel
# - Perl "rename" function

# Citations

# Wu, Y. W., Simmons, B. A. & Singer, S. W. 
# MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. 
# Bioinformatics 32, 605-607, doi:10.1093/bioinformatics/btv638 (2016).

# Kang, D. D., Froula, J., Egan, R. & Wang, Z. 
# MetaBAT, an efficient tool for accurately reconstructing single genomes from complex microbial communities. 
# PeerJ 3, e1165, doi:10.7717/peerj.1165 (2015).

# Anaconda Software Distribution. Computer software. Continuum Analytics. Web. <https://continuum.io>.

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool.

# ------------------------------------------------------------------------------

# Arguments
asm=$1
map_dir=$2
out_dir=$3
threads=$4


mkdir -p $out_dir/maxbin_cov_profiles/ $out_dir/bins_maxbin2/

# Activate virtual environment
source activate binmate

# Note
# MaxBin2 tries to resolve the location of its scripts through a Perl function that does not handle symbolic links properly
# Solution: use the complete absolute path to the MaxBin2 perl script 
# A similar problem is observed with FragGeneSan - Prepend FragGeneScan to PATH to avoid problems
export PATH=/usr/local/bioinfo/FragGeneScan:$PATH

# Generate coverage profiles (MetaBAT script)
jgi_summarize_bam_contig_depths --outputDepth $out_dir/maxbin_cov_profiles/depth_metabat.tsv --minContigDepth 2 $map_dir/*.sorted.bam

index=$(seq 4 2 $(($(ls -l $map_dir/*.sorted.bam | wc -l) * 2 + 2)))
samples=$(ls $map_dir/*sorted.bam | xargs -n 1 basename | sed 's/\..*//')
parallel --link "cut -f 1,{1} $out_dir/maxbin_cov_profiles/depth_metabat.tsv | sed -e '1d' > $out_dir/maxbin_cov_profiles/{2}.tsv" ::: $index ::: $samples

rm $out_dir/maxbin_cov_profiles/depth_metabat.tsv
find $out_dir/maxbin_cov_profiles/ -name "*.tsv" > $out_dir/maxbin_cov_profiles/cov_profiles.txt

# Run MaxBin 2.0
$(readlink -f `which run_MaxBin.pl`) -contig $asm -abund_list $out_dir/maxbin_cov_profiles/cov_profiles.txt -thread $threads -out $out_dir/maxbin2

# Move and rename bins
mv $out_dir/*.fasta $out_dir/bins_maxbin2/
cd $out_dir/bins_maxbin2/
rename "s/fasta/fa/" *.fasta

# Deactivate virtual environment
source deactivate
