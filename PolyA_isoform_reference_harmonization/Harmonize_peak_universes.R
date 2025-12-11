library(dplyr)
library(data.table)
library(tibble)
library(tidyverse)
library(GenomicRanges)
library(stringr)

# directory of separate cell subpopulation references
ref_path <- "path_to_directory_of_cell_subpopulation_peak_universes"
ref_dir <- dir(ref_path)

genes_cat <- data.frame()
for (ref in ref_dir){
  print(paste("Extracting genes from ref:",ref,sep=" "))
  this_ref <- paste(ref_path,ref,sep="\\")
  this_ref <- read.table(this_ref,header=TRUE,sep="\t")
  this_ref_genes <- as.data.frame(this_ref$gene)
  genes_cat <- rbind(genes_cat,this_ref_genes)
}

genes_unique <- unique(genes_cat)
colnames(genes_unique) <- "unique_genes"

ref_cat <- data.frame()
genes_counted <- 0
for (gene in genes_unique$unique_genes){
  gene_cat <- data.frame()
  for (ref in ref_dir) {
    this_ref <- paste(ref_path,ref,sep="\\")
    this_ref <- read.table(this_ref,header=TRUE,sep="\t")
    gene_subset <- this_ref[this_ref$gene==gene,]
    gene_subset <- gene_subset[order(gene_subset$start),]
    gene_cat <- rbind(gene_cat,gene_subset)
  }
  gene_tu <- unique(gene_cat$tu)
  gene_granges <- makeGRangesFromDataFrame(gene_cat,keep.extra.columns=TRUE,ignore.strand=FALSE)
  gene_reduced <- reduce(gene_granges)
  gene_reduced <- as.data.frame(gene_reduced)  
  if (length(gene_reduced$seqnames)>1){
    if (unique(gene_reduced$strand)=="-"){
      peak_cts <- paste("P",length(gene_reduced$start):1,sep="")
      final_annot <- paste(gene_tu,gene,peak_cts,sep=":")   
    }
    else {
      peak_cts <- paste("P",1:length(gene_reduced$start),sep="")
      final_annot <- paste(gene_tu,gene,peak_cts,sep=":")  
    }
  }
  else {
    peak_cts <- "P0"
    final_annot <- paste(gene_tu,gene,peak_cts,sep=":")
  }
  gene_reduced$final_annotation <- final_annot
  gene_reduced$gene <- gene
  ref_cat <- rbind(ref_cat,gene_reduced)
  genes_counted <- genes_counted+1
  if (genes_counted%%100 == 0){print(paste("Completed ",genes_counted," genes",sep=""))}
}

# save the harmonized reference as a txt file
file_out <- paste(getwd(),"Harmonized_Peak_Reference.txt",sep="\\")
write.table(ref_cat,file_out,quote=FALSE,sep="\t",row.names=FALSE)

## END