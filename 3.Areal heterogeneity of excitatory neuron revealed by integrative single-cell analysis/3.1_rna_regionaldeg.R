library(Seurat)
c4_ex <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2011/1109/c4_ex.rds")
DefaultAssay(c4_ex) <- "RNA"

layer <- c("L23it","L45it","L5it","L6it","l56itcar3","l56np","L5et","l6b","l6ct")
for (i in layer) {
  ct <- subset(c4_ex,subset = layer_2 == i)
  ct <- NormalizeData(ct)
  Idents(ct)<- ct$region_2
  ct.markers <- FindAllMarkers(ct, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.01)
  write.table(ct.markers,file= paste0("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2101/0115/deg/",i,".csv"),sep=",")
  
}

coembed <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1217/rna_atac_coembed.rds")
meta <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2101/0116/coembed_meta_atac.rds")
meta <- meta[rownames(coembed@meta.data),]
coembed@meta.data <- meta
c4 <- subset(coembed,subset = batch %in% c("a0723","b0909"))

celltype <- c("LAMP5","PVALB","RELN","SST","VIP","AST","ENDO","PERI","MIC","OLI","OPC")

for (i in celltype) {
  ct <- subset(c4,subset = celltype == i)
  ct <- NormalizeData(ct)
  Idents(ct)<- ct$region_2
  ct.markers <- FindAllMarkers(ct, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.01)
  write.table(ct.markers,file= paste0("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2101/0115/deg/",i,".csv"),sep=",")
}



