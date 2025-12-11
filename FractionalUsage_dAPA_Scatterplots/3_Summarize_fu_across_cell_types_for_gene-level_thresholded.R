library(data.table)
library(dplyr)
library(stringr)
library(Matrix)
library(MASS)

data_dir <- paste0(getwd(),"\\cell_type_fractional_usages\\")
out_dir <- paste0(getwd(),"\\cell_type_fractional_usages\\")

# first assess for median filtered zeros excluded
data_sub_dir_path <- paste0(data_dir,"\\summary_usage_median_filtered_zeros_excluded_Gene-level\\")
data_sub_dir <- list.files(path=data_sub_dir_path,pattern=".csv")

file_ct<-1
for (file in data_sub_dir) {
  
  which_cell <- str_split(file,"_")[[1]][1]
  col_label1 <- paste0(which_cell,"_Normalized_Isoform_Counts")
  col_label2 <- paste0(which_cell,"_Fractional_Usage")
  file_path <- paste0(data_sub_dir_path,file)
  this_data <- read.csv(file_path,row.names=1)
  this_data <- this_data[order(this_data$Isoform_IDs),]
  colnames(this_data)[1] <- col_label1
  colnames(this_data) <- gsub("\\.","-",colnames(this_data))
  
  data_reorder <- this_data[,c(2,3,1,4)]
  which_cols <- c(col_label1,col_label2)
  data_subset <- data_reorder[,which_cols]
  
  if (file_ct==1){data_out <- data_reorder} else {data_out <- cbind(data_out,data_subset)}
  file_ct<-file_ct+1
}

file_out <- paste0(out_dir,"summary_usage_median_filtered_zeros_excluded_Gene-level_All_Cells.csv")
write.csv(data_out,file_out)

## END