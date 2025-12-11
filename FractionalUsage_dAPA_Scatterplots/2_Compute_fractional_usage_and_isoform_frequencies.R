library(data.table)
library(dplyr)
library(stringr)
library(Matrix)
library(MASS)

data_dir <- paste0(getwd(),"\\cell_type_isoform_distributions\\")
out_dir <- paste0(getwd(),"\\cell_type_fractional_usages\\")
if(!dir.exists(out_dir)) {dir.create(out_dir)}

## GENE-LEVEL FILTERED
# first assess for median filtered without zeros 
data_sub_dir_path <- paste0(data_dir,"\\summary_counts_median_filtered_zeros_excluded_Gene-level\\")
data_sub_dir <- list.files(path=data_sub_dir_path,pattern="_zeroed_counts_norm.csv")

data_out_dir <- paste0(out_dir,"\\summary_usage_median_filtered_zeros_excluded_Gene-level\\")
if(!dir.exists(data_out_dir)) {dir.create(data_out_dir)}

for (file in data_sub_dir) {
  
  which_cell <- str_split(file,"_")[[1]][1]
  col_label <- paste0(which_cell,"_Fractional_Usage")
  file_path <- paste0(data_sub_dir_path,file)
  this_data <- read.csv(file_path,row.names=1)
  this_data <- this_data[order(this_data$Isoform_IDs),]
  this_data[ ,col_label]=NA
  
  for (gene in unique(this_data$Genes)) {
    
    this_gene <- this_data[this_data$Genes==gene,]
    total_counts <- sum(this_gene$Normalized_Isoform_Counts)
    if (total_counts==0){
      frac_counts<-replicate(length(this_gene$Isoform_IDs),NA)
      this_data[this_data$Genes==gene,4]=frac_counts}
    else {
      frac_counts <- this_gene$Normalized_Isoform_Counts/total_counts
      this_data[this_data$Genes==gene,4]=frac_counts
    }
  }
  which_summary <- str_split(file,".csv")[[1]][1]
  file_out <- paste0(data_out_dir,which_summary,"_fractional_usage.csv")
  write.csv(this_data,file_out)
}

## END