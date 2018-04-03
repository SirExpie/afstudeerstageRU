#!/bin/bash
set -o errexit;

# Utility: create bin tsv file
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Description
# Generate a tab-seperated file of the format:
# <contig_id>\t<bin_assignment>
# Such a file can be used as input for DAS Tool

# Usage
# util_bin2tsv.sh <bin_dir> <file_extension> <out_file>

# ------------------------------------------------------------------------------

# Arguments
bin_dir=$1
ext=$2
out_file=$3


# Generate file
grep ">" $bin_dir/*$ext* | sed -e 's@.*/@@' | sed -e 's/\.[^.]*\:>/\t/' | awk '{print $2"\t"$1}' > $out_file
