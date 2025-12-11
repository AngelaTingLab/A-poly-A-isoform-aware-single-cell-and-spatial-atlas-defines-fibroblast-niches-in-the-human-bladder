library(Seurat)
library(data.table)
library(stringr)
library(igraph)
library(ggplot2)
library(pheatmap)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(tidyverse)

# load xenium object
xenium_obj <- readRDS("path_to_xenium_object.rds"))

# extract centroids from Xenium object
centroid_data <- GetTissueCoordinates(xenium_obj)

# extract cell types from xenium object metadata, check order with barcodes, append to centroid data
cell_type_data <- as.data.frame(xenium_obj$final_cell_types)
colnames(cell_type_data) <- "cell_types"
cell_type_data$barcodes <- row.names(cell_type_data)
cell_type_data <- cell_type_data[match(centroid_data$cell,cell_type_data$barcodes),]
centroid_data$cell_types <- cell_type_data$cell_types

# initiate neighborhood analysis
nearest5_cells <- list()
n_select <- 5 # how many neighbors to consider
for (i in 1:length(centroid_data$cell_types)){
  
  this_cell <- centroid_data[i,]
  all_other_cells <- centroid_data[-i, ]
  distances <- sqrt((all_other_cells$x - this_cell$x)^2 + (all_other_cells$y - this_cell$y)^2)
  all_other_cells$distances <- distances
  nearest_cells <- all_other_cells[order(all_other_cells$distances), ][1:n_select, ]
  nearest5_cells[[i]] <- nearest_cells
  
  if (i %% 100 == 0){
    print(i)
  }
}

combined_nearest_points <- do.call(rbind,nearest5_cells)
nearest_cell_types <- as.data.frame(combined_nearest_points$cell_types)
nearest_cell_types <- as.data.frame(matrix(nearest_cell_types$`combined_nearest_points$cell_types`,ncol=n_select,byrow=TRUE))
colnames(nearest_cell_types) <- c("Nearest1","Nearest2","Nearest3","Nearest4","Nearest5")
nearest_cell_types$reference_cell <- centroid_data$cell_types

file_out <- paste0("Xenium_NearestNeighbor_Analysis_",n_select,"NN.csv")
write.csv(nearest_cell_types,file_out,row.names=FALSE)

nn_freq <- data.frame()
for (cell in unique(nearest_cell_types$reference_cell)){
  
  cell_subset <- nearest_cell_types[nearest_cell_types$reference_cell==cell,]
  nn_vector <- as.data.frame(cell_subset[1:n_select])
  nn_vector <- unlist(nn_vector)
  nn_vector <- data.frame(value=nn_vector)
  nn_frequency <- as.data.frame(table(nn_vector$value))
  nn_frequency$reference_cell <- cell
  colnames(nn_frequency)[1] <- "neighbors"

  nn_frequency <- nn_frequency %>%
    mutate(Percent=(Freq/sum(Freq))*100,y_position=cumsum(Percent)-0.5*Percent) %>%
    ungroup()
  nn_freq <- rbind(nn_freq,nn_frequency)
}

file_out <- paste0(out_dir,"Xenium_NearestNeighbor_Analysis_",n_select,"NN_Frequency_Table.csv")
write.csv(nn_freq,file_out,row.names=FALSE)

## END

