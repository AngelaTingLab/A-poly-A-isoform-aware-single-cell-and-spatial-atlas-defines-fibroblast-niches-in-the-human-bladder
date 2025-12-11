library(data.table)
library(dplyr)
library(stringr)
library(Matrix)
library(ggplot2)
library(viridis)
library(patchwork)
library(MASS)
library(DescTools)
library(gghighlight)

# functions
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# script starts
data_dir <- paste0(getwd(),"\\cell_type_fractional_usages\\")
out_dir <- paste0(getwd(),"\\cell_type_dAPA\\")
if(!dir.exists(out_dir)) {dir.create(out_dir)}

# from now on, we are only using the median filtered gene-level threshold isoform dataset
# this is the dataset that was selected after optimization of thresholding parameters
umap_colors <- c("cell_type_0"="#006666","cell_type_1"="#FFC679","cell_type_2"="#66CCCC","cell_type_3"="#990099","cell_type_4"="#FF00FF",
                 "cell_type_5"="#FF99FF","cell_type_6"="#3333FF","cell_type_7"="#0099FF","cell_type_8"="#9933FF")

plot_out_dir <- paste0(out_dir,"cell_type_dAPA_plots\\")
file_out_dir <- paste0(out_dir,"cell_type_dAPA_gene_candidate_lists\\")
if(!dir.exists(plot_out_dir)) {dir.create(plot_out_dir)}
if(!dir.exists(file_out_dir)) {dir.create(file_out_dir)}

data_path <- paste0(data_dir,"\\summary_usage_median_filtered_zeros_excluded_Gene-level_All_Cells.csv")
data <- read.csv(data_path,row.names=1)
colnames(data) <- gsub("\\.","-",colnames(data))
gene_type <- as.data.frame(str_split_fixed(data$Isoform_IDs,":",3))
data <- data[gene_type$V3 %in% c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11"),]

frac_usage_only <- data[,seq(3,ncol(data),2)]
cell_types <- as.data.frame(str_split_fixed(colnames(frac_usage_only),"_",4))
cell_types <- as.data.frame(cell_types$V1)
colnames(cell_types) <- "cell_types"
cell_pairs <- as.data.frame(t(as.data.frame(combn(cell_types$cell_types,2,simplify=FALSE))))
row.names(cell_pairs)<-NULL

for (pair_num in (1:length(cell_pairs$V1))){
  
  cell1 <- cell_pairs$V1[pair_num]
  cell2 <- cell_pairs$V2[pair_num]
  
  cell1_count_label <- paste0(cell1,"_Normalized_Isoform_Counts")
  cell1_frac_label <- paste0(cell1,"_Fractional_Usage")
  
  cell2_count_label <- paste0(cell2,"_Normalized_Isoform_Counts")
  cell2_frac_label <- paste0(cell2,"_Fractional_Usage")
  
  cell1_counts <- data[,c("Isoform_IDs","Genes",cell1_count_label)]
  cell1_frac_usage <- data[,c("Isoform_IDs","Genes",cell1_frac_label)]
  
  cell2_counts <- data[,c("Isoform_IDs","Genes",cell2_count_label)]
  cell2_frac_usage <- data[,c("Isoform_IDs","Genes",cell2_frac_label)]
  
  cell1_has <- subset(cell1_frac_usage,(!is.na(cell1_frac_usage[,3])))
  cell2_has <- subset(cell2_frac_usage,(!is.na(cell2_frac_usage[,3])))
  both_have <- intersect(cell1_has$Isoform_IDs,cell2_has$Isoform_IDs)
  
  cell1_frac_use <- cell1_has[cell1_has$Isoform_IDs %in% both_have,]
  cell2_frac_use <- cell2_has[cell2_has$Isoform_IDs %in% both_have,]
  
  cell1_count_use <- cell1_counts[cell1_counts$Isoform_IDs %in% both_have,]
  cell2_count_use <- cell2_counts[cell2_counts$Isoform_IDs %in% both_have,]
  
  data_to_plot <- cbind(cell1_count_use,cell1_frac_use[,3],cell2_count_use[,3],cell2_frac_use[,3])
  colnames(data_to_plot)[4]<-cell1_frac_label
  colnames(data_to_plot)[5]<-cell2_count_label
  colnames(data_to_plot)[6]<-cell2_frac_label
  data_to_plot$cell_frac_dist <- abs((data_to_plot[,6]-data_to_plot[,4]))
  data_to_plot$frac_fc <- abs(log2(data_to_plot[,6]/data_to_plot[,4]))
  data_to_plot$APA_Key <- NA
  apa_cell1_label <- paste0(cell1," Enriched")
  apa_cell2_label <- paste0(cell2," Enriched")
  thresh <- quantile(data_to_plot$cell_frac_dist,probs=c(0.95))
  thresh1 <- 0.15
  thresh2 <- log2(1.5)
  
  data_to_plot$APA_Key[((data_to_plot$cell_frac_dist)>=thresh1)&(data_to_plot$frac_fc>=thresh2)&(data_to_plot[,4]>data_to_plot[,6])]=apa_cell1_label
  data_to_plot$APA_Key[((data_to_plot$cell_frac_dist)>=thresh1)&(data_to_plot$frac_fc>=thresh2)&(data_to_plot[,4]<data_to_plot[,6])]=apa_cell2_label
  data_to_plot$APA_Key[((data_to_plot$cell_frac_dist)<thresh1)|((data_to_plot$frac_fc)<thresh2)]="Comparable Expression"
  data_to_plot$APA_Key <- factor(data_to_plot$APA_Key,levels=c("Comparable Expression",apa_cell1_label,apa_cell2_label))
  
  cell1_color <- unname(umap_colors[cell1])
  cell2_color <- unname(umap_colors[cell2])
    
  apa_scatter <- ggplot(data_to_plot,aes(x=data_to_plot[,4],y=data_to_plot[,6],color=APA_Key)) + geom_point() + theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="bottom") + xlab(cell1_frac_label) + ylab(cell2_frac_label) +
    scale_color_manual(values=c("#CCCCCC",cell1_color,cell2_color))
  apa_scatter

  # save the apa scatter
  plot_out <- paste0(plot_out_dir,cell1,"_v_",cell2,"_APA.jpg")
  ggsave(plot_out,apa_scatter,width=20,height=20,units="cm",dpi=600)
  
  # save the isoform ids for extreme apa candidates
  candidates_to_save <- data_to_plot[data_to_plot$APA_Key %in% c(apa_cell1_label,apa_cell2_label),]
  candidates_to_save$Mean_Normalized_Counts <- (candidates_to_save[,3]+candidates_to_save[,5])/2
  candidates_to_save <- candidates_to_save[order(-candidates_to_save$Mean_Normalized_Counts),]
  file_out <- paste0(file_out_dir,cell1,"_v_",cell2,"_APA_candidates.csv")
  write.csv(candidates_to_save,file_out)
  
  # save all normalized count, fractional usage, and average counts per comparison
  data_to_plot$Mean_Normalized_Counts <- (data_to_plot[,3]+data_to_plot[,5])/2
  col_label <- paste0(cell1,"_v_",cell2,"_Mean_Normalized_Counts")
  colnames(data_to_plot)[10] <- col_label
  data_to_plot$cell_frac_dist <- NULL
  data_to_plot$frac_fc <- NULL
  data_to_plot$APA_Key <- NULL
  file_out <- paste0(file_out_dir,cell1,"_v_",cell2,"_ALL_isoforms.csv")
  write.csv(data_to_plot,file_out)
  
}

## END