library(Seurat)
library(data.table)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(ReactomePA)

# load deg list
deg_list <- read.csv("path_to_deg_file.csv"))
gene_entrez_ids <- mapIds(org.Hs.eg.db,keys=fb_dapa_genes[,1],column="ENTREZID",keytype="SYMBOL",multiVals="first")

# KEGG terms, dotplot, and saving
kegg_enrichment <- enrichKEGG(gene=gene_entrez_ids,organism='hsa',pvalueCutoff=0.05)
kegg_summary <- as.data.frame(kegg_enrichment)
kegg_dotplot <- dotplot(kegg_enrichment,showCategory=20)

write.csv(kegg_summary,paste0(out_dir,"\\KEGG_Terms_",which_fb,".csv"))
ggsave(paste0(out_dir,"\\KEGG_Terms_DotPlot_",which_fb,".jpg"),plot=kegg_dotplot,dpi=600,width=8,height=10)

# run analysis with go
# BP terms, dotplot, and saving
go_enrichment <- enrichGO(gene=gene_entrez_ids,OrgDb=org.Hs.eg.db,ont="BP",pvalueCutoff=0.05)
go_summary <- as.data.frame(go_enrichment)
go_dotplot <- dotplot(go_enrichment,showCategory=20)
write.csv(go_summary,paste0(out_dir,"\\GO_BP_Terms_",which_fb,".csv"))
ggsave(paste0(out_dir,"\\GO_BP_Terms_DotPlot_",which_fb,".jpg"),plot=go_dotplot,dpi=600,width=8,height=10)

# CC terms, dotplot, and saving
go_enrichment <- enrichGO(gene=gene_entrez_ids,OrgDb=org.Hs.eg.db,ont="CC",pvalueCutoff=0.05)
go_summary <- as.data.frame(go_enrichment)
go_dotplot <- dotplot(go_enrichment,showCategory=20)
write.csv(go_summary,paste0(out_dir,"\\GO_CC_Terms_",which_fb,".csv"))
ggsave(paste0(out_dir,"\\GO_CC_Terms_DotPlot_",which_fb,".jpg"),plot=go_dotplot,dpi=600,width=8,height=10)

# MF terms, dotplot, and saving
go_enrichment <- enrichGO(gene=gene_entrez_ids,OrgDb=org.Hs.eg.db,ont="MF",pvalueCutoff=0.05)
go_summary <- as.data.frame(go_enrichment)
go_dotplot <- dotplot(go_enrichment,showCategory=20)
write.csv(go_summary,paste0(out_dir,"\\GO_MF_Terms_",which_fb,".csv"))
ggsave(paste0(out_dir,"\\GO_MF_Terms_DotPlot_",which_fb,".jpg"),plot=go_dotplot,dpi=600,width=8,height=10)

# ReactomePA terms, dotplot, and saving
rpa_enrichment <- enrichPathway(gene=gene_entrez_ids,organism="human",pvalueCutoff=0.05)
rpa_summary <- as.data.frame(rpa_enrichment)
rpa_dotplot <- dotplot(rpa_enrichment,showCategory=20)
write.csv(rpa_summary,paste0(out_dir,"\\ReactomePA_Terms_",which_fb,".csv"))
ggsave(paste0(out_dir,"\\ReactomePA_Terms_DotPlot_",which_fb,".jpg"),plot=rpa_dotplot,dpi=600,width=8,height=10)

## END

