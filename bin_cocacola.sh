#!/bin/bash
set -o errexit;

# Metagenome binning - COCACOLA
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# bin_cocacola.sh <contigs.fasta> <map_dir> <out_dir> <threads>

# Dependencies
# - COCACOLA
# - CONCOCT
# - Anaconda with "binmate" virtual environment
# - perl "rename" function

# Citations

# Lu, Y. Y., Chen, T., Fuhrman, J. A. & Sun, F. 
# COCACOLA: binning metagenomic contigs using sequence COmposition, read CoverAge, CO-alignment and paired-end read LinkAge. 
# Bioinformatics 33, 791-798, doi:10.1093/bioinformatics/btw290 (2017).

# Alneberg, J. et al. 
# Binning metagenomic contigs by coverage and composition. 
# Nature methods 11, 1144-1146, doi:10.1038/nmeth.3103 (2014).

# Anaconda Software Distribution. Computer software. Continuum Analytics. Web. <https://continuum.io>.

# ------------------------------------------------------------------------------

# Arguments
asm=$1
map_dir=$2
out_dir=$3
threads=$4

# Activate virtual environment
source activate binmate

mkdir -p $out_dir/bins_cocacola

# Generate coverage profiles (CONCOCT script)
bam_files=$(find $map_dir -name "*.sorted.bam")
gen_input_table.py $asm $bam_files | cut -f1,3- > $out_dir/cov_cocacola.tsv

# Generate composition data (CONCOCT script)
fasta_to_features.py $asm $(grep -c ">" $asm) 4 $out_dir/kmer_cocacola.tsv

# Generate linkage data (CONCOCT script)
find $map_dir -name "*.sorted.bam" -exec basename {} \; > $out_dir/samples.txt
bam_to_linkage.py -m $threads --regionlength 500 --fullsearch --samplenames $out_dir/samples.txt $asm $bam_files > $out_dir/linkage_concoct.tsv
$(dirname $(readlink -f $0))/util_edges_concoct2cocacola.py -i $out_dir/linkage_concoct.tsv -o $out_dir/edges_cocacola.tsv

# Run COCACOLA
# NOTE: The cocacola.py script has been modified.
# The original script calls to other scripts in its "auxiliary" folder using os.getcwd() (present working dir)
# By modifying line 196 - 199: "os.getcwd()" to "os.path.dirname(os.path.realpath(__file__))",
# it is no longer required to have a copy of COCACOLA in the working directory.
python $(which cocacola.py) --contig_file $asm --abundance_profiles $out_dir/cov_cocacola.tsv --composition_profiles $out_dir/kmer_cocacola.tsv --edge_list $out_dir/edges_cocacola.tsv --output $out_dir/bins_cocacola.txt

# Extract clusters/bins (CONCOCT script)
extract_fasta_bins.py --output_path $out_dir/bins_cocacola/ $asm $out_dir/bins_cocacola.txt

# Move and rename bins
cd $out_dir/bins_cocacola/
rename 's/^/cocacola\./' *.fa

# Deactivate virtual environment
source deactivate
