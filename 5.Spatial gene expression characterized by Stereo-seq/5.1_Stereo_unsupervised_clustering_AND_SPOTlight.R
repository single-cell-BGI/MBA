library(Seurat)
library(dplyr)
library(SPOTlight)
library(RColorBrewer)
library(Matrix)

### unsupervised clustering
data_bin1 = data.table::fread("cortex_M1_bin1.txt", data.table = F)

data_matrix$coor_x = round(data_bin1$x / 50)
data_matrix$coor_y = round(data_bin1$y / 50)

data_matrix$bin = paste0(data_matrix$coor_x, "_", data_matrix$coor_y)

data_matrix2 = data_matrix %>% group_by(bin, geneID) %>% summarise(UMIs = sum(MIDCounts))
counts = sparseMatrix(as.numeric(as.factor(data_matrix2$geneID)), as.numeric(as.factor(data_matrix2$bin)),x=data_matrix2$UMIs)
rownames(counts) = levels(as.factor(data_matrix2$geneID))
colnames(counts) = levels(as.factor(data_matrix2$bin))

M1_SP <- CreateSeuratObject(counts = counts, project = "Spatial", min.cells = 10, min.features = 10)
M1_SP[["percent.mt"]] <- PercentageFeatureSet(M1_SP, pattern = "^mt-")

M1_SP@meta.data$coor_x <-  unlist(lapply(strsplit(split="_", rownames(M1_SP@meta.data)),'[[',1))
M1_SP@meta.data$coor_y <-  unlist(lapply(strsplit(split="_", rownames(M1_SP@meta.data)),'[[',2))

M1_SP <- SCTransform(M1_SP, verbose = FALSE)

M1_SP <- RunPCA(M1_SP, verbose = FALSE, assay = 'SCT')
M1_SP <- RunUMAP(M1_SP, dims = 1:15, verbose = FALSE)
M1_SP <- FindNeighbors(M1_SP, dims = 1:15, verbose = FALSE)
M1_SP <- FindClusters(M1_SP, verbose = FALSE)

saveRDS(M1_SP, file = "Cortex_M1_SP.RDS")

###
###
### Processed the single nucleus data

cortex_sc_M1 <- readRDS("Cortex_M1.RDS")

set.seed(123)
cortex_sc_M1 <- SCTransform(cortex_sc_M1, verbose = FALSE) %>% RunPCA(., verbose = FALSE) %>% RunUMAP(., dims = 1:30, verbose = FALSE)

Idents(cortex_sc_M1) <- cortex_sc_M1$CellType

cluster_markers_all <- FindAllMarkers(object = cortex_sc_M1, 
                                              min.pct = 0.4, 
                                              verbose = FALSE, 
                                              only.pos = TRUE)

###################################################################################
### SPOTlight: deconvolution the cell type in each bin #############################

set.seed(123)

spotlight_ls <- spotlight_deconvolution(se_sc = cortex_sc_M1,
                                        counts_spatial = M1_SP@assays$RNA@counts,
                                        clust_vr = CellType,
                                        cluster_markers = cluster_markers_all,
                                        cl_n = 50, # 100 by default
                                        hvg = 3000,
                                        ntop = NULL,
                                        transf = "uv",
                                        method = "nsNMF",
                                        min_cont = 0.09)


decon_mtrx <- spotlight_ls[[2]]
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

decon_mtrx <- decon_mtrx[, which(colnames(decon_mtrx) != "res_ss")]

decon_mtrx[decon_mtrx < 0.2 ] = 0

M1_SP@meta.data <- cbind(M1_SP@meta.data, decon_mtrx)

predict_CellType = apply(decon_mtrx, 1, function(x){
                                   index = which.max(x)
                                   a = colnames(decon_mtrx)[index]
                                   return(a)})

M1_SP$CellType = predict_CellType

saveRDS(M1_SP, file = "Cortex_M1_SP_SPOTlight.rds")

#################################################
#### DEG of cell type in the Stereo-seq #########

Idents(M1_SP) <- M1_SP$predict_CellType

cluster_markers_sp <- FindAllMarkers(object = M1_SP, 
                                              min.diff.pct = 0.2, 
                                              verbose = FALSE, 
                                              only.pos = TRUE)

write.csv(cluster_markers_sp, file = "Cortex_M1_SP_CellType_DEG.csv", quote = F)
