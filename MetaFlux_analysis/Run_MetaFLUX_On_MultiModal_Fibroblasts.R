library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(METAFlux)
library(stringi)
library(stringr)
library(tidyr)
library(osqp)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#output directory for metabolite flux and other data/plots
out_dir <- paste0(getwd(),"\\METAFlux_Outputs\\")
dir.create(out_dir,showWarnings=FALSE)

# custom functions
generate_boots <- function(celltype, n) {
  dt <- data.frame(cluster = celltype, id = 1:length(celltype))
  index <- do.call(cbind, sapply(1:n, function(x) {
    splits <- dt %>%
      group_by(cluster) %>%
      sample_n(dplyr::n(), replace = TRUE) %>%
      ungroup() %>%
      dplyr::select("id")
  }))
  return(index)
}

get_ave_exp <- function(i, myseurat, samples,myident) {
  meta.data=myseurat@meta.data[samples[,i],]
  sample <-GetAssayData(myseurat, assay = "RNA")[,samples[,i]]
  name <- colnames(sample)
  for (j in 1:length(name)) {
    name[j] <- paste0(name[j], "_", j)
  }
  colnames(sample) <- name
  rownames(meta.data) <- name
  SeuratObject<-suppressWarnings(
    CreateSeuratObject(count=sample, meta.data = meta.data))
  SeuratObject<-NormalizeData(SeuratObject,verbose = TRUE)
  ave<-GetAssayData(AverageExpression(SeuratObject,group.by = myident,return.seurat = T), assay = "RNA") %>% as.data.frame()
  return(ave)
}

edit_calculate_avg_exp <- function(myseurat,myident,n_bootstrap,seed) {
  set.seed(seed)
  samples=generate_boots(myseurat@meta.data[,myident],n_bootstrap)
  exp <- lapply(1:n_bootstrap,get_ave_exp,myseurat,samples,myident)
  exp <- do.call(cbind, exp)
  return(exp)
}

compute_sc_flux <- function(num_cell, fraction, fluxscore, medium) {
  if (sum(fraction) != 1)
    stop("Sum of fractions must be eqaul to 1")
  if (length(fraction) != num_cell)
    stop("Number of cell clusters does not match with the length of fraction")
  # perform preparation
  Hgem <- METAFlux:::Hgem
  mat <- Hgem$S
  reaction_name <- Hgem$Reaction
  names(reaction_name) <- NULL
  A_combined <- METAFlux:::A_combined
  D <- as.sparse(matrix(0, 8378, 13082))
  candi_list <- lapply(rapply(list(mat, lapply(1:2, function(x) D)),
                              enquote, how = "unlist"), eval)
  matrix_construct <- diag(1, ncol = num_cell, nrow = num_cell)
  matrix_construct[matrix_construct == 0] <- 2
  celltype_matrix <- NULL
  A_matrix <- NULL
  for (i in 1:num_cell) {
    A_matrix[[i]] <- A_combined
    celltype_matrix[[i]] <- do.call(cbind, candi_list[matrix_construct[,
                                                                       i]])
  }
  message("Preparing for TME S matrix.....")
  final_s <- rbind(do.call(cbind, A_matrix), do.call(rbind, celltype_matrix))
  whole3 <- rbind(as.sparse(diag(-1, 1648, 1648)), as.sparse(matrix(0,
                                                                    num_cell * nrow(mat), 1648)))
  final_s <- cbind(final_s, whole3)
  external <- paste("external_medium", reaction_name[which(Hgem$pathway ==
                                                             "Exchange/demand reactions")])
  reaction_name[which(Hgem$pathway == "Exchange/demand reactions")] <- paste0("internal_medium ",
                                                                              reaction_name[which(Hgem$pathway == "Exchange/demand reactions")])
  cell_reaction <- NULL
  for (i in 1:num_cell) {
    cell_reaction[[i]] <- paste(paste("celltype", i), reaction_name)
  }
  construct_reaction_names <- c(unlist(cell_reaction), external)
  P1 <- as.sparse(diag(1, ncol(final_s), ncol(final_s)))
  message("S matrix completed......")
  flux_vector <- list()
  message("Compute metabolic flux......")
  Seq <- seq(1, ncol(fluxscore), by = num_cell)
  pb <- txtProgressBar(0, length(Seq), style = 3)
  for (t in Seq) {
    setTxtProgressBar(pb, match(t, Seq))
    score <- fluxscore[, c(t:(t + num_cell - 1))]
    P <- as.sparse(diag(c(unlist(Map(rep, fraction, 13082)), rep(1,
                                                                 1648)), ncol(final_s), ncol(final_s)))
    fraction_finals <- final_s %*% P
    q <- rep(0, ncol(final_s))
    q[seq(13015, ncol(final_s), by = 13082)] <- -10000 * fraction
    A <- as.sparse(rbind(fraction_finals, P1))
    ras <- c(as.vector(unlist(score[, 1:num_cell])), rep(1, 1648))
    origlb <- c(rep(Hgem$LB, num_cell), rep(-1, 1648))
    origlb[rep(Hgem$rev == 1, num_cell)] <- (-ras[rep(Hgem$rev == 1,
                                                      num_cell)])
    origlb[rep(Hgem$rev == 0, num_cell)] <- 0
    origub <- ras
    origlb[tail(1:ncol(final_s), 1648)] <- 0
    matches <- sapply(medium$ID, function(x) intersect(which(stri_detect_fixed(construct_reaction_names,x)), tail(1:ncol(final_s), 1648)))
    matches <- unlist(matches)
    origlb[matches] <- -1
    l <- c(rep(0, nrow(final_s)), origlb)
    u <- c(rep(0, nrow(final_s)), origub)
    settings <- osqpSettings(max_iter = 1000000L, eps_abs = 1e-04,
                             eps_rel = 1e-04, verbose = FALSE, adaptive_rho_interval = 50)
    model <- osqp(P, q, A, l, u, settings)
    # Solve problem
    res <- model$Solve()
    flux_vector <- append(flux_vector, list(res$x))
  }
  close(pb)
  flux_vector <- as.data.frame(do.call(cbind, flux_vector))
  rownames(flux_vector) <- construct_reaction_names
  return(flux_vector)
  
}

# load in multimodal fibroblast data 
cell_colors <- c("sn-3"="#990099","sn-1"="#FFC679","sn-0"="#006666","sn-2"="#66CCCC","sn-4"="#FF00FF","sn-6"="#3333FF","sn-8"="#9933FF","sn-7"="#0099FF","sn-5"="#FF99FF")

seurat_obj <- readRDS("path_to_seurat_object_with_multimodal_fibroblast_labels.rds"))
Idents(seurat_obj) <- seurat_obj$cell_types_multimodal
seurat_obj@active.assay <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
multimodal_umap <- DimPlot(seurat_obj)
multimodal_umap

# run MetaFLUX analysis
#check cell types, here we want to model cell-type level flux
table(seurat_obj$cell_types_multimodal)

#Calculate the mean expression for bootstrap samples
mean_exp <- edit_calculate_avg_exp(myseurat=seurat_obj,myident='cell_types_multimodal',n_bootstrap=50,seed=5)

#calculate metabolic reaction scores
scores<-calculate_reaction_score(data=mean_exp)

#calculate the fractions of celltype/clusters
round(table(seurat_obj$cell_types_multimodal)/nrow(seurat_obj@meta.data),3)
sum(table(seurat_obj$cell_types_multimodal)/nrow(seurat_obj@meta.data)) # should sum to 1
cell_fractions <- unlist(list(table(seurat_obj$cell_types_multimodal)/nrow(seurat_obj@meta.data)))

# num_cell=number of cell types/clusters, here we have 9 cell types, thus input is 9. 
# make sure the sum of cell percentage is 1.
# the order of fraction should be the same as that of "Mean_exp" data
data("human_gem")
flux <- compute_sc_flux(num_cell=9,fraction=c(0.24,0.18,0.16,0.12,0.11,0.09,0.05,0.03,0.02),fluxscore=scores,medium=human_gem)
file_out <- paste0(out_dir,"\\Flux.csv")
write.csv(flux,file_out)

# optional: flux scores normalization
cbrt <- function(x) {sign(x) * abs(x)^(1/3)}
flux=cbrt(flux)
dim(flux)
file_out <- paste0(out_dir,"\\Normalized_Flux.csv")
write.csv(flux,file_out)

# calculate mean flux and make a heatmap
# compute mean flux across bootstraps and reshape for heatmap generation
mean_flux <- as.data.frame(rowMeans(flux))
pathway<-unique(unlist(human_gem$SUBSYSTEM))
cell_types <- c("celltype 1","celltype 2","celltype 3","celltype 4","celltype 5","celltype 6","celltype 7","celltype 8","celltype 9")
all_cells_pathway_scores <- matrix()
for (cell in cell_types){
  print(cell)
  cell_flux_subset <- flux[grepl(cell, rownames(flux)), ]
  pathway_score<-list()
  for (i in pathway){
    path=i
    activity_score<-c()
    for (d in 1:ncol(flux)){
      activity_score[d]<-mean(abs(cell_flux_subset[which(unlist(human_gem$SUBSYSTEM)==i),d]))
    } 
    pathway_score[[i]]<-activity_score
  }
  all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))
  all_pathway_score <- as.data.frame(rowMeans(all_pathway_score))
  colnames(all_pathway_score) <- cell
  all_cells_pathway_scores <- cbind(all_cells_pathway_scores,all_pathway_score)
}

all_cells_pathway_scores$all_cells_pathway_scores <- NULL
colnames(all_cells_pathway_scores) <- c("sn-0","sn-1","sn-2","sn-3","sn-4","sn-5","sn-6","sn-7","sn-8")
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pathway_heatmap <- pheatmap(all_cells_pathway_scores,cluster_cols=TRUE,color = rev(mapal),scale = "row")

file_out <- paste0(out_dir,"\\Pathway_Scores_Heatmap_Mean.jpg")
ggsave(file_out,plot=pathway_heatmap,width=20,height=30,dpi=600)

## END
