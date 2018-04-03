#!/bin/bash
set -o errexit;

# Subsampling of merged quality trimmed paired-end libraries
# ------------------------------------------------------------------------------

# Author: P. Boersma - P.Boersma@science.ru.nl
# Date: November 9, 2017

# Used for
# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Description
# This script will run seqkit command sample which subsamples a (gzipped) fastq file
# for a given number or proportion.

# Usage
# seqkit sample [-p|--proportion] float [-s|--rand-seed] int <fastq file> -o <outname>

# Dependencies
# - seqkit
# - GNU parallel

# Citations
# W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q
# file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962.

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool
# ------------------------------------------------------------------------------

# Arguments
props=$1
out_dir=$2
fq_dir=$3
merge=$4
subsampling=$5
data_dir=$6

proplist="${props//","/" "}"

p=0.10

echo "Merge: "$merge" Sub: "$subsampling

# Merge quality-trimmed paired-end libraries
if [ "$merge" = true ] ; then
  echo "Merging..."
  parallel "cat $fq_dir/*_[rR]{}* > $out_dir/merged_R{}.fq.gz" ::: 1 2
  # subsampling merged dataset
  if [ "$subsampling" = true ] ; then
#    echo "Subsampling with "$proplist"..."
#    for p in $proplist
#    do
      #mkdir $out_dir/sub_$p/
      echo "subsampling R1..."
      seqkit sample -p $p --rand-seed 11 $out_dir/merged_R1.fq.gz -o $out_dir/sub_$p\_merged_R1.fq.gz
      echo "subsampling R2..."
      seqkit sample -p $p --rand-seed 11 $out_dir/merged_R2.fq.gz -o $out_dir/sub_$p\_merged_R2.fq.gz
#    done
  fi

# Subsampling all samples separately
elif [ "$subsampling" = true ] ; then
  echo "Subsampling with "$proplist"..."
  echo $data_dir
  file_list=($(ls $data_dir/*.fq.gz))
  echo ${file_list[@]}
  for file in ${file_list[@]}
  do
    echo $(basename $file)
    for p in $proplist
      do
      if [ ! -d $out_dir\sub_$p/ ]; then
        echo "Making sub-directory"
        mkdir $out_dir\sub_$p/
      fi
      seqkit sample -p $p --rand-seed 11 $file -o $out_dir/sub_$p/sub_$p\_$(basename $file)
      done
  done
fi

# Generate sample stats after sub-sampling
echo "Generating stats..."
seqkit stats $out_dir/*.fq.gz -a | tee $out_dir/asm_stats_seqkit

rm -f $out_dir/merged_R*.fq.gz
