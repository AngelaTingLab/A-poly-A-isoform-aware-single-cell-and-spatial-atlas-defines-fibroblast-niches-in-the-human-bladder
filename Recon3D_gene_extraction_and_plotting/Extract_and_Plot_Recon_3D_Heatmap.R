library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(tidyverse)
library(stringi)
library(stringr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# create output directory
out_dir <- paste0(getwd(),"\\Recon_3D_DB_Analysis_and_Heatmap\\")
dir.create(out_dir,showWarnings=FALSE)

# load in multimodal fibroblast data 
cell_colors <- c("sn-3"="#990099","sn-1"="#FFC679","sn-0"="#006666","sn-2"="#66CCCC","sn-4"="#FF00FF","sn-6"="#3333FF","sn-8"="#9933FF","sn-7"="#0099FF","sn-5"="#FF99FF")

seurat_obj <- readRDS("path_to_multimodal_fibroblast_object.rds")
multimodal_umap <- DimPlot(seurat_obj)
seurat_obj$cell_types_multimodal <- factor(seurat_obj$cell_types_multimodal,
                                           levels=c("sn-4","sn-3","sn-1","sn-2","sn-0","sn-6","sn-7","sn-5","sn-8"))
Idents(seurat_obj) <- seurat_obj$cell_types_multimodal

# set assay to RNA, normalize, and scale data
seurat_obj@active.assay <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# load the recon 3D DB
recon_3d <- read.csv("path_to_Recon_3D_DB_Extracted.csv")
recon_3d_genes <- as.data.frame(unique(recon_3d$Gene_Name))
recon_3d_genes <- as.data.frame(recon_3d_genes[-1,])
colnames(recon_3d_genes) <- "Gene"

# subset for the Recon 3D genes
available_genes <- intersect(recon_3d_genes$Gene, rownames(seurat_obj))
subset_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[available_genes, ]
normalized_df <- as.data.frame(subset_matrix)
normalized_df <- as.data.frame(t(normalized_df))
normalized_df$cell_types_multimodal <- seurat_obj$cell_types_multimodal

# run DEA and add pct1-pct2 for these genes
recon_obj <- CreateSeuratObject(counts = t(normalized_df[,1:1818]))
metadata_df <- data.frame(cell_types_multimodal = normalized_df$cell_types_multimodal)
metadata_df$cell_ids <- row.names(normalized_df)
rownames(metadata_df) <- colnames(recon_obj) 
recon_obj <- AddMetaData(recon_obj, metadata = metadata_df)
Idents(recon_obj) <- recon_obj$cell_types_multimodal

recon_obj <- NormalizeData(recon_obj)
recon_obj <- ScaleData(recon_obj)
mm_fb_recon_markers <- FindAllMarkers(recon_obj)
mm_fb_recon_markers$'pct1-pct2' <- (mm_fb_recon_markers$pct.1)-(mm_fb_recon_markers$pct.2)
mm_fb_recon_markers <- mm_fb_recon_markers[mm_fb_recon_markers$p_val_adj<0.05,]
mm_fb_recon_markers <- mm_fb_recon_markers[mm_fb_recon_markers$avg_log2FC>log2(1.5),]

file_out <- paste0(out_dir,"\\Recon_3D_DEGs_No_Pct_Filt.csv")
write.csv(mm_fb_recon_markers,file_out)

mm_fb_recon_markers_pct_delt_01 <- mm_fb_recon_markers[mm_fb_recon_markers$`pct1-pct2`>=0.1,]
file_out <- paste0(out_dir,"\\MM_FB_Recon_3D_DEGs_0.1_Pct_Filt.csv")
write.csv(mm_fb_recon_markers_pct_delt_01,file_out)

mm_fb_recon_markers_pct_delt_005 <- mm_fb_recon_markers[mm_fb_recon_markers$`pct1-pct2`>=0.05,]
file_out <- paste0(out_dir,"\\MM_FB_Recon_3D_DEGs_0.05_Pct_Diff_Filt.csv")
write.csv(mm_fb_recon_markers_pct_delt_005,file_out)

# compute the mean and median expression for each gene
mean_expr_by_celltype <- normalized_df %>%
  group_by(cell_types_multimodal) %>%
  summarise(across(all_of(available_genes), mean, na.rm = TRUE))

# subset for cell types that you want to plot, here we plot the layered fibroblasts
mean_expr_by_celltype <- mean_expr_by_celltype[colnames(mean_expr_by_celltype) %in% c("sn-0","sn-1","sn-2","sn-3","sn-4","sn-6")]
mean_expr_by_celltype <- mean_expr_by_celltype[,c("sn-4","sn-3","sn-1","sn-2","sn-0","sn-6")]

# compute and optionally implement low variance filter
var_thresh <- 0.001
data_var <- apply(mean_expr_by_celltype,1,var)
data_filt <- mean_expr_by_celltype[data_var > var_thresh, , drop = FALSE]
data_filt <- as.matrix(data_filt)
storage.mode(data_filt) <- "numeric"
data_to_plot <- log2(data_filt+1)

# plot a heatmap of these genes and cluster it
cell_colors <- c("sn-4"="#FF00FF","sn-3"="#990099","sn-1"="#FFC679","sn-2"="#66CCCC",
                 "sn-0"="#006666","sn-6"="#3333FF")
col_annotation <- data.frame(Cell_Type = colnames(data_filt))
rownames(col_annotation) <- colnames(data_filt)
col_annotation$Cell_Type <- factor(col_annotation$Cell_Type,
                                   levels=c("sn-4","sn-3","sn-1","sn-2","sn-0","sn-6"))
cell_type_colors <- setNames(cell_colors, col_annotation$Cell_Type)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
recon_3D_heatmap <- pheatmap(data_to_plot,cluster_cols=FALSE,cluster_rows=TRUE,fontsize_row=4,
                             color=rev(mapal),scale="row",clustering_distance_rows="euclidean",
                             border_color=NA,annotation_col=col_annotation,annotation_colors=list(Cell_Type=cell_type_colors))

file_out <- paste0(out_dir,"\\Recon_3D_Heatmap_LP_FB_Only_Var_filtered_Var",var_thresh,".jpg")
ggsave(file_out,plot=recon_3D_heatmap,width=10,height=10,dpi=600)

scaled_mat <- t(scale(t(data_to_plot)))
rescaled_matrix <- scaled_mat[row_order,]
file_out <- paste0(out_dir,"\\Recon_3D_Heatmap_LP_FB_Only_Var_filtered_Var",var_thresh,"_RESCALED_MATRIX.csv")
write.csv(rescaled_matrix,file_out)

row_order <- recon_3D_heatmap$tree_row$order
col_order <- recon_3D_heatmap$tree_col$order
clustered_matrix <- data_to_plot[row_order, col_order]

file_out <- paste0(out_dir,"\\Recon_3D_Heatmap_LP_FB_Only_Var_filtered_Var",var_thresh,".csv")
write.csv(clustered_matrix,file_out)

## END
