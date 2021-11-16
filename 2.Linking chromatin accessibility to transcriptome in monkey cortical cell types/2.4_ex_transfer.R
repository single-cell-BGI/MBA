library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)

coembed <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1217/rna_atac_coembed.rds")
rna <- subset(coembed,subset = celltype_4 %in% c("EX1","EX2","EX3","EX4","EX5","EX6","EX7","EX8","EX9","EX10","EX11","EX12","EX13","EX14") & batch2 == "RNA")
atac <- subset(coembed,subset = celltype_4 %in% c("EX1","EX2","EX3","EX4","EX5","EX6","EX7","EX8","EX9","EX10","EX11","EX12","EX13","EX14") & batch2 == "atac")

DefaultAssay(atac) <- "ACTIVITY"
atac <- FindVariableFeatures(atac)
atac <- NormalizeData(atac)
atac <- ScaleData(atac)

DefaultAssay(atac) <- "ATAC"
VariableFeatures(atac) <- names(which(Matrix::rowSums(atac) > 100))
atac <- RunLSI(atac, n = 50, scale.max = NULL)
atac <- RunUMAP(atac, reduction = "lsi", dims = 1:50)

rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna, resolution = 0.5)
rna <- RunUMAP(rna, dims = 1:10)
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
layer.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$layer_2, 
                                  weight.reduction = atac[["lsi"]])
atac <- AddMetaData(atac, metadata = layer.predictions)
saveRDS(atac@meta.data,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1223/atac_ex_layer_meta")
atac$layer_2 = atac$predicted.id
saveRDS(atac@meta.data,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1223/atac_ex_layer_meta")

saveRDS(atac,file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1223/atac_ex_layer.rds")
pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1223/atac_layer.pdf",width = 12,height = 5)
p1 <- DimPlot(atac, group.by = "layer_2", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(rna, group.by = "layer_2", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
plot_grid(p1,p2)
dev.off()


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

coembed$layer_2 <- ifelse(!is.na(coembed$layer_2), coembed$layer_2, coembed$predicted.id)


saveRDS(coembed, file = "/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1223/ex_layer.rds")

pdf("/hwfssz5/ST_PRECISION/TOMCAT/UCB_share/syn/macabrain/2012/1223/ex_layer.pdf")
p1 <- DimPlot(coembed, reduction = "umap", group.by = "layer_2",label = T,repel = T)
print(p1)
dev.off()
