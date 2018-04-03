#!/bin/bash
set -o errexit;

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl


# SETUP
# ------------------------------------------------------------------------------

# Load configuration file
dos2unix -q $1
source $1

# Activate virtual environment
source activate binmate

# Resolve the location of BinMate
SOURCE="${BASH_SOURCE[0]}"

while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
binmate="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# Project root dir
proj_dir=$(pwd)/binmate_$prefix
echo "Project directory: "$proj_dir

# Create project folder structure
mkdir -p $proj_dir/{dat,qc/{pre_trim,post_trim},trim,asm/{asm_stats,ctg_stats},map/{se,pe},tax/contigs,bin/{binsanity,cocacola,concoct,maxbin2,metabat2,das_tool},\
checkm/{binsanity,cocacola,concoct,maxbin2,metabat2,das_tool},results}
echo "Folders created..."

# Symlink data to project folder
find $fq_dir -maxdepth 1 -type f -regextype posix-extended -iregex '.*/.*\.fastq.*|.*/.*\.fq.*' -exec ln -s {} $proj_dir/dat/ \;

cd $proj_dir

# Save config
cat $1 >> .conf_$$

# SEQUENCE DATA QUALITY PROCESSING
# ------------------------------------------------------------------------------

# FASTQC (pre-trimming)
#echo "Generating FASTQC report pre-trimming..."
#$binmate/qc_fastqc.sh dat/ qc/pre_trim/

# Quality based trimming of paired-end reads - BBTools BBDuk
echo "Trimming..."
$binmate/trim_pe_bbduk.sh dat/ $ref $k_bbduk $mink $hdist $ktrim $qtrim $trimq $maq $ftl $ftr $maxns $minlen trim/ $threads

# FASTQC (post-trimming)
#echo "Generating FASTQC report post-trimming..."
#$binmate/qc_fastqc.sh trim/ qc/post_trim/


# SEQUENCE CLASSIFICATION
# ------------------------------------------------------------------------------

# Classify reads - Kaiju
if [ "$classify_reads" = true ] ; then
  echo "Classify reads..."
  mkdir -p tax/reads/
  parallel "cat trim/*_[rR]{}* > trim/merged_R{}.fq.gz" ::: 1 2
  $binmate/classify_pe_reads_kaiju.sh all_pe_reads $kaiju_db trim/merged_R1.fq.gz trim/merged_R2.fq.gz tax/reads/ $threads
  rm trim/merged*
fi

# METAGENOME COASSEMBLY
# ------------------------------------------------------------------------------

# Merging and subsizing fastq files
if [ "$merge" = true ] || [ "$subsampling" = true ] ; then
  echo "Merging and subsizing..."
  $binmate/sub_sampling.sh $subsize asm/ trim/ $merge $subsampling dat/
fi
echo "=====================EAT SLEEP PROGRAM REPEAT============================"
sleep 60

# De novo assembly
if [ "$metaspades" = true ] ; then
  $binmate/asm_metaspades.sh $spades_k $spades_correct asm/ $memory $threads $merge $subsampling
elif [ "$idba_ud" = true ] ; then
  $binmate/asm_idba_ud.sh trim/ asm/ $threads
elif [ "$megahit" = true ] ; then
  $binmate/asm_megahit.sh trim/ asm/ $threads
fi

# Generate assembly and contig statistics
$binmate/qc_asm_stats.py -i asm/*/*_contigs.fasta --report --datatable -o asm/

# Filter on contig size
$binmate/qc_asm_stats.py -i asm/*/*_contigs.fasta --report --datatable --length $min_ctg_len --fasta -o asm/


# READ MAPPING
# ------------------------------------------------------------------------------

asm=asm/*_contigs*_gt$min_ctg_len.fasta

# Map PE + SE reads to metagenome - BWA
$binmate/map_bwa.sh $asm trim/ map/ $threads
rm map/*unsorted.bam map/*.sam
mv map/*singles* map/se/
mv map/*bam* map/pe/


# METAGENOME BINNING
# ------------------------------------------------------------------------------

# BinSanity
if [ $(grep -c ">" $asm) -lt 75000 ] ; then
  cd bin/binsanity/
  $binmate/bin_binsanity.sh ../../$asm ../../map/pe/ $threads
  cd $proj_dir
else
  rm -r bin/binsanity/ checkm/binsanity/
fi

# COCACOLA
$binmate/bin_cocacola.sh $asm map/pe/ bin/cocacola/ $threads

# CONCOCT
$binmate/bin_concoct.sh $asm map/pe/ bin/concoct/ $threads

# MaxBin 2.0
$binmate/bin_maxbin2.sh $asm map/pe/ bin/maxbin2/ $threads

# MetaBAT
$binmate/bin_metabat2.sh $asm $min_ctg_len map/pe/ bin/metabat2/ $threads


# CONSENSUS BINNING
# ------------------------------------------------------------------------------

# Generate bin file input for DAS Tool
bin_files=bin/das_tool/bin_files/
mkdir -p $bin_files

parallel "$binmate/util_bin2tsv.sh {} fa $bin_files/{/}.tsv" ::: $(find bin/*/ -maxdepth 1 -type d -name "bins_*")

# Generate argument strings
das_tool_files=$(find $bin_files -name "*.tsv" -printf %p, | sed -e 's/.$//')
das_tool_labels=$(echo $das_tool_files | sed -e "s:${bin_files}bins_::g;s:\.tsv::g")

# Run DAS TOOL
$binmate/bin_das_tool.sh $asm $das_tool_files $das_tool_labels bin/das_tool/ $threads


# BIN QUALITY CONTROL
# ------------------------------------------------------------------------------

# CheckM
bin_labels=$(find bin/ -maxdepth 1 -mindepth 1 -type d | sed -e 's/.*\///' | sort)
bin_dirs=$(find bin/*/ -maxdepth 1 -type d -name "bins_*" | sort)
checkm_out_dir=$(find checkm/* -maxdepth 1 -type d | sort)

parallel -j 4 --link "$binmate/qc_checkm.sh {1} $asm {2} fa map/pe/ {3} 5" ::: $bin_labels ::: $bin_dirs ::: $checkm_out_dir


# COMPILE FINAL RESULTS
# ------------------------------------------------------------------------------

# Gather results
ln -s asm/contigs.fasta results/                                # Contigs assembly
ln -s $asm results/                                             # Contigs assembly (filtered)
ln -s bin/das_tool/bins_das_tool/ results/                      # FASTA files final bins (DAS Tool)
ln -s checkm/das_tool/das_tool.bin_qa_plot.pdf results/         # CheckM plot (DAS Tool)
ln -s checkm/das_tool/das_tool_16s_seq/ results/                # 16S rRNA sequences (DAS Tool)
ln -s checkm/das_tool/das_tool.unbinned_contigs.fasta results/  # FASTA file unbinned contigs (DAS Tool)

# Generate binning summary tables + plots
checkm_results=$(find checkm/*/ -type f -name "checkm_results*" -printf %p, | sed -e 's/.$//')
labels=$(find checkm/*/ -type f -name "checkm_results*" -printf %f, | sed -e 's/checkm_results\.//g;s/\.tsv//g;s/.$//')
Rscript $binmate/qc_bin_eval.R $checkm_results $labels results/

# Generate contig information table + bin plots
# $binmate/util_ctg_table.sh $asm map/pe/ bin/das_tool/bins_das_tool/ results/ $threads


# Deactivate virtual environment
source deactivate

