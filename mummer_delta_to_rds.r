### Scripts covert Delta files from mummer alignments adapted from:
# https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/
# https://github.com/Uauy-Lab/pangenome-haplotypes/blob/master/haplotype-block-assignment/assign_mummer_blocks_whole_genome.r

library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)


readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

calculate_perc_id_mid_points <- function(data){
  data$r_length <- (data$re - data$rs)
  data$perc_id <- ((data$r_length - data$error)/data$r_length)*100
  data$perc_id_factor <- data$perc_id
  data[data$perc_id < 100, "perc_id_factor"] <- "<100"
  data$perc_id_factor <- as.factor(data$perc_id_factor)
  data$r_mid <- (data$rs + data$re)/2
  data$q_mid <- (data$qs + data$qe)/2
  return(data)
}

pre_plot_analysis <- function(delta_path, min_size = 20000){
  data = readDelta(delta_path)
  data <- calculate_perc_id_mid_points(data)
  data_filt <- data[data$r_length >= min_size,]
  return(data_filt)
}

base_dir <- "/Volumes/quirozj/09_watseq/15_alignments_to_ibspy/whole_genome_mummer/"
varieties_to_plot <- list.files(paste0(base_dir, "aln"))
chromosomes_to_plot <- list("chr1A", "chr2A", "chr3A", "chr4A", "chr5A", "chr6A","chr7A")
print(chromosomes_to_plot)

#we will only use the 20kb filter here (this is already filtered using mummer)
min_size <- 20000

for (chrom in chromosomes_to_plot){
  chrom_aln <- data.frame("rs" = numeric(),
                      "re" = numeric(),
                      "qs" = numeric(),
                      "qe" = numeric(),
                      "error" = numeric(),
                      "qid" = character(),            
                      "rid" = character(),            
                      "strand" = character(),         
                      "r_length" = numeric(),
                      "perc_id" = numeric(),
                      "perc_id_factor" = factor(),
                      "r_mid" = numeric(),
                      "q_mid" = numeric(),
                      "comparison" = character(),
                      "chrom" = character())

  for (variety in  varieties_to_plot){
    variety_dir <- paste0(base_dir, "aln/", variety)
    chrom_dir <- paste0(variety_dir, "/", chrom)
    all_files <- list.files(chrom_dir)

    #This gets just the filtered deltas 
    filtered_delta <- all_files[grep("filtered_L20Kb_rq.delta", all_files)]
    for (comparison in filtered_delta){
      print(comparison)
      comparison_delta_path <- paste0(chrom_dir, "/", comparison)
      
      comparison_filt <- pre_plot_analysis(delta_path = comparison_delta_path, min_size = min_size)
      comparison_filt$comparison <- comparison
      comparison_filt$chrom <- chrom
      chrom_aln <- rbind(chrom_aln, comparison_filt)
    }
  }
  saveRDS(chrom_aln, file = paste0("/Volumes/quirozj/09_watseq/15_alignments_to_ibspy/mummer_rds_mon/all_20_kb_filtered_delta_", chrom, "_tables.rds"))
}