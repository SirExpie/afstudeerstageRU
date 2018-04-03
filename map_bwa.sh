#!/bin/bash
set -o errexit;

# Map reads to metagenome - BWA
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# map_bwa.sh <contigs.fasta> <trim_fq_dir> <out_dir> <threads>
# Expects paired-end reads (trimmed) in de-interleaved format (forward and reverse reads contained in separate fastq files).
# Files may be compressed (gzip). Prefixes of forward and reverse file should be identical (e.g. sample1_R1.fq.gz, sample1_R2.fq.gz).

# Dependencies
# - BWA Burrows-Wheeler Aligner
# - SAMtools (version 1.3.1 using htslib 1.3.1)
# - GNU parallel

# Citations

# Li, H. & Durbin, R. 
# Fast and accurate long-read alignment with Burrows-Wheeler transform. 
# Bioinformatics 26, 589-595, doi:10.1093/bioinformatics/btp698 (2010).

# Li, H. et al. 
# The Sequence Alignment/Map format and SAMtools. 
# Bioinformatics 25, 2078-2079, doi:10.1093/bioinformatics/btp352 (2009).

# O. Tange (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
asm=$1
fq_dir=$2
out_dir=$3


# Determine number of threads per job
threads=$(expr $4 / $(ls $fq_dir/*R1* | wc -l))

# Function: map PE reads to a reference using BWA
function map_bwa_pe {
  asm=$1
  r1=$2
  r2=$3
  out_dir=$4
  threads=$5
  echo $asm $r1 $r2 $out_dir $threads

  prefix=$(echo $2 | sed -e 's:.*/::;s:_[rR][12].*$::')
  base=$(basename $asm)
  echo "Prefix: "$prefix
  echo "Base: "$base
  bwa mem -M -t $threads $asm $r1 $r2 > $out_dir/$base\_$prefix.sam
  samtools fixmate -@ $threads -O bam $out_dir/$base\_$prefix.sam $out_dir/$base\_$prefix.unsorted.bam
  samtools sort -o $out_dir/$base\_$prefix.sorted.bam -O bam -@ $threads $out_dir/$base\_$prefix.unsorted.bam 
  samtools index -b -@ $threads $out_dir/$base\_$prefix.sorted.bam
}

# Function: map SE reads to a reference using BWA
function map_bwa_se {
asm=$1
s1=$2
out_dir=$3
threads=$4

prefix=$(echo $2 | sed -e 's:.*/::;s:\..*$::')
base=$(basename $asm)
echo "Prefix: "$prefix
echo "Base: "$base
bwa mem -M -t $threads $asm $s1 > $out_dir/$base\_$prefix.sam
samtools view -@ $threads -bt $asm -o $out_dir/$base\_$prefix.unsorted.bam $out_dir/$base\_$prefix.sam
samtools sort -@ $threads $out_dir/$base\_$prefix.unsorted.bam -o $out_dir/$base\_$prefix.sorted.bam
samtools index -b -@ $threads $out_dir/$base\_$prefix.sorted.bam
}

export -f map_bwa_pe
export -f map_bwa_se

# Index assembly
bwa index $asm

# Run PE mapping jobs
parallel --link "map_bwa_pe $asm {1} {2} $out_dir $threads" ::: $(find -L $fq_dir -maxdepth 1 -type f -regextype posix-extended -iregex '.*/.*_R1.*' | sort) ::: $(find -L $fq_dir -maxdepth 1 -type f -regextype posix-extended -iregex '.*/.*_R2.*' | sort)

# Run SE mapping jobs
parallel --link "map_bwa_se $asm {1} $out_dir $threads" ::: $(find -L $fq_dir -maxdepth 1 -type f -regextype posix-extended -iregex '.*/.*_singles_.*' | sort)
