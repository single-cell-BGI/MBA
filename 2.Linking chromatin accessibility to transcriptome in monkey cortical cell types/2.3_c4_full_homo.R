library(Seurat)
library(ggplot2)
library(cowplot)


all <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1010/homo_maca_ex_inte_layer.rds")

DefaultAssay(all) <- "RNA"

c4 <-subset(all,subset = batch %in% c("a0723","b0909"))
full <- subset(all,subset = batch == "full_length")
homo <- subset(all,subset = batch == "homo")

c4_homo.list <- list(c4,homo)

c4_homo.list <- lapply(X = c4_homo.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

c4_homo.anchors <- FindIntegrationAnchors(object.list = c4_homo.list, dims = 1:20)
c4_homo.combined <- IntegrateData(anchorset = c4_homo.anchors, dims = 1:20)
c4_homo.combined <- ScaleData(c4_homo.combined, verbose = FALSE)
c4_homo.combined <- RunPCA(c4_homo.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
c4_homo.combined <- RunUMAP(c4_homo.combined, reduction = "pca", dims = 1:20)
c4_homo.combined <- FindNeighbors(c4_homo.combined, reduction = "pca", dims = 1:20)
c4_homo.combined <- FindClusters(c4_homo.combined, resolution = 0.5)

saveRDS(c4_homo.combined,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1012/c4_homo_inte.rds")

pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1012/c4_homo_cluster.pdf",height = 5,width = 12)
p1 <- DimPlot(c4_homo.combined, reduction = "umap", group.by = "species")
p2 <- DimPlot(c4_homo.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
DefaultAssay(c4_homo.combined) <- "RNA"

pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1012/c4_homo_species.pdf",height = 26, width =26)
p3 <- DimPlot(c4_homo.combined, reduction = "umap", split.by = "species")
print(p3)
dev.off()

full_homo.list <- list(full,homo)

full_homo.list <- lapply(X = full_homo.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

full_homo.anchors <- FindIntegrationAnchors(object.list = full_homo.list, dims = 1:20)
full_homo.combined <- IntegrateData(anchorset = full_homo.anchors, dims = 1:20)
full_homo.combined <- ScaleData(full_homo.combined, verbose = FALSE)
full_homo.combined <- RunPCA(full_homo.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
full_homo.combined <- RunUMAP(full_homo.combined, reduction = "pca", dims = 1:20)
full_homo.combined <- FindNeighbors(full_homo.combined, reduction = "pca", dims = 1:20)
full_homo.combined <- FindClusters(full_homo.combined, resolution = 0.5)

saveRDS(full_homo.combined,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1012/full_homo_inte.rds")

pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1012/full_homo_cluster.pdf",height = 5,width = 12)
p1 <- DimPlot(full_homo.combined, reduction = "umap", group.by = "species")
p2 <- DimPlot(full_homo.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
DefaultAssay(full_homo.combined) <- "RNA"

pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1012/full_homo_species.pdf",height = 26, width =26)
p3 <- DimPlot(full_homo.combined, reduction = "umap", split.by = "species")
print(p3)
dev.off()


