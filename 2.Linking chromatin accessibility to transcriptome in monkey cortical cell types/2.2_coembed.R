library(ArchR)
library(Seurat)
library(SummarizedExperiment)
library(cowplot)
library(ggplot2)
atac <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1206/atac_umap.rds")
#1206跑的umap过后的atac矩阵
rna <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2010/1014/all_0917_celltype.rds")
meta <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1216/all_meta_color.rds")
ameta <- meta[rownames(rna@meta.data),]
rna@meta.data <- meta
DefaultAssay(rna) <- "RNA"
rna$batch2 = "RNA"
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna,resolution = 0.5)
rna <- RunUMAP(rna, dims = 1:10)
saveRDS(rna,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1206/rna_umap.rds")

transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$celltype_4, 
                                     weight.reduction = atac[["lsi"]])
atac <- AddMetaData(atac, metadata = celltype.predictions)
table(atac$prediction.score.max > 0.5)
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]])

# this line adds the imputed data matrix to the atac object
atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- FindNeighbors(coembed, dims = 1:10)
coembed <- FindClusters(coembed, resolution = 0.5)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype_4 <- ifelse(!is.na(coembed$celltype_4), coembed$celltype_4, coembed$predicted.id)
saveRDS(coembed@meta.data,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1217/rna_atac_coembedmeta.rds")
saveRDS(coembed,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1217/rna_atac_coembed.rds")
pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1217/coembed_celltype.pdf",height = 5,width = 14)
p1 <- DimPlot(coembed, group.by = "batch2")
p2 <- DimPlot(coembed, group.by = "celltype_4", label = TRUE, repel = TRUE)
plot_grid(p1,p2)
dev.off()
pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1217/coembed_cluster.pdf")
p1 <- DimPlot(coembed,reduction = "umap", label = TRUE)
print(p1)

marker<-c("SLC17A7","GAD1","RELN","VIP","LAMP5","PVALB","SST","MOG","PDGFRA","PDGFRB","APBB1IP","FLT1","SLC1A2")
pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1217/coembed_marker.pdf",height = 20,width = 20)
p <- FeaturePlot(coembed, features = marker, order = T)
print(p)
dev.off()


