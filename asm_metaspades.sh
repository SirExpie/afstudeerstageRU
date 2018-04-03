#!/bin/bash
set -o errexit;

# Co-assembly of metagenome - metaSPAdes
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Description
# - This script will run metaSPAdes using all available paired-end libraries by merging them together.
# - Single-end libraries generated during trimming of PE reads are not used.
# - Error correction prior to assembly is recommended, however, this process requires a lot of memory (~10 Gb for 1 E. coli genome).
# - It may be necessary to skip the error correction step to prevent SPAdes from running out of memory and aborting.
# - If the assembler runs out of memory, consider normalizng read coverage (bbnorm, BBTools)

# Usage
# asm_metaspades.sh <trim_fq_dir> <kmers> correct [true|false] <out_dir> <mem(gb)> <threads>

# Dependencies
# - MetaSPAdes
# - BBtools (optional)
# - GNU parallel

# Citations

# Nurk, S., Meleshko, D., Korobeynikov, A. & Pevzner, P. A.
# metaSPAdes: a new versatile metagenomic assembler.
# Genome research, doi:10.1101/gr.213959.116 (2017).
# ------------------------------------------------------------------------------

# Arguments
kmers=$1
correct=$2
asm_dir=$3
mem=$4
threads=$5
merge=$6
subsamp=$7

# De novo metagenome assembly
if [ "$correct" = true ]; then
  # Error correction + assembly
  echo "Assembling..."
  if [ "$merge" = true ] || [ "$subsamp" = true ]; then
    echo $(ls -f $asm_dir/*.fq.gz)
    for file in $(ls -f $asm_dir/*.fq.gz) ;
    do
      echo "Left reads file: "$file
      if [[ $file == *"_R1"* ]]
      then
        file2=$(echo $file | sed 's/_R1/_R2/g')
        echo "Right reads file: "$file2
        outdir=$(echo ${file%.fq.gz})
        echo "Outdir: "$outdir
        spades.py --meta -k $kmers --threads $threads --memory $mem -o $outdir/ --pe1-1 $file --pe1-2 $file2

        file_out=$(echo ${file%.fq.gz})
        file_out="${file_out}_contigs.fasta"
        echo "Output name: "$file_out
        mv $outdir/contigs.fasta $outdir/$(basename $file_out)
        echo "file and file2: "$file" "$file2
        rm $file || true
        rm $file2 || true
      else
        continue
      fi
    done
  else
    echo $(ls -f trim/*.fq.gz)
    for file in $(ls -f trim/*.fq.gz) ;
    do
      echo "Left reads file: "$file
      if [[ $file == *"_R1"* ]]
      then
        file2=$(echo $file | sed 's/_R1/_R2/g')
        echo "Right reads file: "$file2
        outdir=$(echo ${file%.fq.gz})
        outdir=$(basename $outdir)
        outdir=$(echo "asm/"$outdir)
        echo "Outdir: "$outdir
        spades.py --meta -k $kmers --threads $threads --memory $mem -o $outdir/ --pe1-1 $file --pe1-2 $file2

        file_out=$(echo ${file%.fq.gz})
        file_out="${file_out}_contigs.fasta"
        echo "Output name: "$file_out
        mv $outdir/contigs.fasta $outdir/$(basename $file_out)
      fi
    done
  fi
else
  # Skip error correction
    echo "Assembling..."
  echo $(ls -f $asm_dir*/*.fq.gz)
  for file in $(ls -f $asm_dir/*/*.fq.gz) ;
  do
    echo "Left reads file: "$file
    if [[ $file == *"_R1"* ]]
    then
      file2=$(echo $file | sed 's/_R1/_R2/g')
      echo "Right reads file: "$file2
      outdir=$(echo ${file%.fq.gz})
      echo "Outdir: "$outdir
      spades.py --only-assembler --meta -k $kmers --threads $threads --memory $mem -o $outdir/ --pe1-1 $file --pe1-2 $file2

      file_out=$(echo ${file%.fq.gz})
      file_out="${file_out}_contigs.fasta"
      echo "Output name: "$file_out
      mv $outdir/contigs.fasta $asm_dir/$(basename $file_out)
      rm -f $file || true
      rm -f $file2 || true
    else
      continue
    fi
  done
fi
