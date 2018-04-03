# Reformat CheckM tables and generate plots
# ------------------------------------------------------------------------------

# BinMate - Metagenome binning pipeline
# J .Frank - j.frank@science.ru.nl

# Usage
# qc_eval_bin.R <checkm_table_i.tsv,*table_n> <label_i,label_n> <out_dir>

# Table should contain the joined results (use "checkm join_tables") of:
# - checkm tree_qa (phylogeny markers)
# - checkm qa (bin quality, should also include coverage information generated using "checkm coverage")
# - checkm profile (mapped reads)

# Dependencies
# - ggplot2

# Citations

# Wickham, H. 
# ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
# ------------------------------------------------------------------------------

library(ggplot2)

args = commandArgs(TRUE)

# Functions

# Load CheckM data, format table
create_table = function(prefix, checkm_file) {
  
  # Read CheckM table
  header = unlist(strsplit(readLines(checkm_file, n=1), split="\t"))
  data = read.table(checkm_file, sep="\t", header=F, stringsAsFactors=F, skip=1)
  
  # Select columns with coverage data
  # Calculate total coverage
  cov_grep = paste("Coverage \\(", collapse="|")
  idx_cov = grep(cov_grep, header, perl=T)[order(grep(cov_grep, header, perl=T, value=T))]
  cov = data[,idx_cov]
  tot_cov = apply(cov, 1, sum)
  
  # Select columns with mapping data
  # Calculate total number of mapped reads
  map_reads_grep = paste("\\: mapped reads", collapse="|")
  idx_map_reads = grep(map_reads_grep, header, perl=T)[order(grep(map_reads_grep, header, perl=T, value=T))]
  idx_map_reads = grep(map_reads_grep, header, perl=T)
  tot_map_reads = apply(data[,idx_map_reads], 1, sum)
  
  perc_reads_grep = paste("\\% mapped", collapse="|")
  idx_perc_reads = grep(perc_reads_grep, header, perl=T)[order(grep(perc_reads_grep, header, perl=T, value=T))]
  idx_perc_reads = grep(perc_reads_grep, header, perl=T)
  perc_reads = data[,idx_perc_reads]
  
  # Select all columns of interest and generate master table
  sel = grep(paste(c("Bin Id", "contained", "sister", "Completeness", "Contamination", "heterogeneity", "Genome size \\(bp\\)", "GC$", "\\# contigs", "N50 \\(contigs\\)"), collapse="|"), header, perl=T)
  data = cbind(prefix, data[,sel], tot_cov, cov, tot_map_reads, perc_reads)
  
  # Generate formatted header
  header = c("method", header[sel], "tot_coverage", header[idx_cov], "tot_mapped_reads", header[idx_perc_reads])
  header = gsub(" ", "_", header)
  header = gsub("[():]", "", header)
  header = gsub("#", "no", header)
  header = gsub("%", "perc", header)
  header = casefold(header)
  colnames(data) = header
  
  data = data.frame(data)
  data = data[order(data$completeness, decreasing=T),] # sort
  return(data)
}

# Output table
write_table = function(data, out_dir){
  # Write table to file
  out_file = sprintf("%s/%s_bin_results.tsv", out_dir, casefold(names(data)))
  write.table(data.frame(data[[1]]), out_file, quote=F, row.names=F, col.names=T, sep="\t")
  return(NULL)
}

create_plot = function(data, type, max_cont, min_comp, out_dir) {
  
  p = ggplot(data, aes(x=contamination, y=completeness)) +
    geom_point(aes(size=tot_coverage, color=n50_contigs), shape=1) +
    scale_color_gradient(low="black", high="#f40ed6", trans="log") +
    labs(x="\nContamination", y="Completeness\n", size="Coverage\n", color="N50") + 
    lims(x=c(0, max_cont), y=c(min_comp, 100)) +
    facet_wrap( ~facet_label) +
    theme_bw() +
    theme(aspect.ratio=1)
  
  out_file = sprintf("%s/p_bining_performance_%s-bins.png", out_dir, type)
  ggsave(filename=out_file, plot=p, width=10, height=10)
}

# ------------------------------------------------------------------------------

checkm_files = unlist(strsplit(args[1], ","))
prefixes = unlist(strsplit(args[2], ","))
out_dir=args[3]
dat_list = list() # Data container

# Quit script if number of files and prefixes are unequal
if(length(checkm_files) != length(prefixes)){
  print("Number of files and prefixes should be the same. Exit.")
  stop()
}

# Parse CheckM files
for(i in seq(1,length(checkm_files))){
  dat_list[[prefixes[i]]] = create_table(prefixes[i], checkm_files[i])
}

# Create table
for(i in seq(1, length(dat_list))) {
  write_table(dat_list[i], out_dir)
}

# Create plot
for(i in seq(1, length(dat_list))) {
  dat_list[[i]]$facet_label = sprintf("%s (n=%s)", dat_list[[i]]$method, nrow(dat_list[[i]]))
}

all_dat = data.frame(do.call("rbind", dat_list))
all_dat$facet_label = factor(all_dat$facet_label, levels=unique(all_dat$facet_label))

create_plot(all_dat, "all", max(all_dat$contamination), 0, out_dir)
create_plot(all_dat, "normal", 50, 0, out_dir)
create_plot(all_dat, "hq", 10, 75, out_dir)
