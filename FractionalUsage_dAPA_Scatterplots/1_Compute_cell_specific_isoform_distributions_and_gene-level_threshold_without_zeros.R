library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(Matrix)
library(MASS)

which_cell <- "cell_type"

out_dir <- paste0(getwd(),"\\cell_type_isoform_distributions\\")
if(!dir.exists(out_dir)) {dir.create(out_dir)}

counts_norm <- read.csv("path_to_normalized_counts.csv")
isoform_ids <- row.names(counts_norm)
cell_ids <- colnames(counts_norm)
print(length(unique(isoform_ids)))

iso_sum_norm <- as.data.frame(rowSums(counts_norm))
iso_sum_norm["Isoform_IDs"] <- rownames(iso_sum_norm)
colnames(iso_sum_norm) <- c("Normalized_Isoform_Counts","Isoform_IDs")
iso_sum_norm <- as.data.frame(iso_sum_norm[order(-iso_sum_norm$Normalized_Isoform_Counts),])
colnames(iso_sum_norm) <- c("Normalized_Isoform_Counts","Isoform_IDs")

sum_norm_genes <- as.data.frame(str_split_fixed(iso_sum_norm$Isoform_IDs,":",3))
iso_sum_norm['Genes'] <- sum_norm_genes$V2
gene_sum_norm <- ddply(iso_sum_norm,"Genes", numcolwise(sum))
gene_sum_norm_no_zeros <- gene_sum_norm[gene_sum_norm$Normalized_Isoform_Counts>0,]

gene_norm_summary <- summary(gene_sum_norm_no_zeros)
poiss_test <- fitdistr(gene_sum_norm_no_zeros$Normalized_Isoform_Counts,"Poisson")
lambda <- poiss_test$estimate
se <- poiss_test$sd

# filter data by median isoform expression and save it
median_filtered_dir <- paste0(out_dir,"\\summary_counts_median_filtered_zeros_excluded_Gene-level")
if(!dir.exists(median_filtered_dir)) {dir.create(median_filtered_dir)}

q2 <-quantile(gene_sum_norm_no_zeros$Normalized_Isoform_Counts,probs=c(0.5))
q2_gene_filt <- gene_sum_norm_no_zeros[gene_sum_norm_no_zeros$Normalized_Isoform_Counts>=q2,]
q2_gene_filt_keep <- unique(q2_gene_filt$Genes)
q2_gene_filt_omit <- setdiff(iso_sum_norm$Genes,q2_gene_filt_keep)

q2_filt_isoforms_filt <- iso_sum_norm[iso_sum_norm$Genes %in% q2_gene_filt_keep,]
q2_filt_isoforms_zeroed <- iso_sum_norm
q2_filt_isoforms_zeroed[q2_filt_isoforms_zeroed$Genes %in% q2_gene_filt_omit,1]=0

q2_out <- paste0(median_filtered_dir,"\\",which_cell,"_gene-level_median_filtered_counts_norm.csv")
write.csv(q2_filt_isoforms_filt,q2_out)

q2_out <- paste0(median_filtered_dir,"\\",which_cell,"_gene-level_median_zeroed_counts_norm.csv")
write.csv(q2_filt_isoforms_zeroed,q2_out)

q2_isoforms <- length(q2_filt_isoforms_filt$Normalized_Isoform_Counts)
q2_genes <- length(unique(q2_filt_isoforms_filt$Genes))
  
print("-------")
print(paste0("Number of Q2 filtered isoforms = ",q2_isoforms))
print("-------")
print(paste0("Number of Q2 filtered genes = ",q2_genes))
print("-------")

## END