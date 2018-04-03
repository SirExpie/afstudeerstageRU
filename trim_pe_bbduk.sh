#!/bin/bash
set -o errexit;

# Quality based trimming of paired-end reads - BBTools: BBDuk
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# trim_pe_bbduk.sh <fq_dir> <ref> <k> <mink> <hdist> <ktrim> <qtrim> <trimq> <maq> <ftl> <ftr> <maxns> <minlen> <out_dir> <threads>

# Dependencies
# - BBTools
# - GNU parallel

# Citations

# Bushnell B. - BBTools 
# bbushnell@lbl.gov
# https://sourceforge.net/projects/bbmap/

# Tange O. (2011): GNU Parallel - The Command-Line Power Tool

# ------------------------------------------------------------------------------

# Arguments
fq_dir=$1

# Calculate number of threads per job
threads=$(expr ${15} / $(ls $fq_dir/*_[rR]1* | wc -l))

# Adapter trimming
ref=$2                          # FASTA file containing sequences to screen against; e.g. adapters, primers, PhiX (spike-ins)
k_bbduk=$3                      # Kmer length used for identifying contaminants (includes reverse-complement searching, contaminants < K cannot be found)
mink=$4                         # Minimum lenght of Kmer used when approaching the end of the read
hdist=$5                        # Number of mismatches allowed in screening for contamination (Hamming distance)
ktrim=$6                        # Enable adapter/contamination trimming; either right-trimming -3' adapters (r), left-trimming -5' adapters (l)
# tbo                           # Flag: trim adapters based on pair overlap detection

# Quality trimming/filtering
qtrim=$7                        # Enable quality trimming (Phred algorithm): either right-trimming -3' (r), left-trimming -5' (l),  both (rl) or sliding window (w)
trimq=$8                        # Average quality trimming treshold
maq=$9                          # Discards reads with an average quality < maq

# Length trimming/filtering
ftl=${10}                       # Trim leftmost bases (positions: 0-tfl)
ftr=${11}                       # Trim the rightmost number of bases (positions: ftr-end)
maxns=${12}                     # Discard reads containing > maxns ambiguous bases
minlen=${13}                    # Discard reads < minlen (bp. post-processing)

out_dir=${14}

# Performing trimming
parallel --link "bbduk.sh in={1} in2={2} \
    ref=$ref k=$k_bbduk mink=$mink hdist=$hdist ktrim=$ktrim tbo qtrim=$qtrim trimq=$trimq maq=$maq ftl=$ftl ftr=$ftr maxns=$maxns minlen=$minlen tossjunk=t t=$threads \
    stats=$out_dir/{=1 s:.*/::;s:_[rR]1::;s:\..*$:_contamination.txt:; =} \
    out=$out_dir/{=1 s:.*/::;s:\..*$:_trim.fq.gz:; =} \
    out2=$out_dir/{=2 s:.*/::;s:\..*$:_trim.fq.gz:; =} \
    outs=$out_dir/{=1 s:.*/::;s:_[rR]1::;s:\..*$:_singles_trim.fq.gz:; =}" ::: $(find -L $fq_dir -maxdepth 1 -type f -regextype posix-extended -iregex '.*/.*_R1.*' | sort) ::: $(find -L $fq_dir -maxdepth 1 -type f -regextype posix-extended -iregex '.*/.*_R2.*' | sort)
