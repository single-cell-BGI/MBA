args=commandArgs(T)
sub <- args[1]

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
.libPaths("/hwfssz5/ST_PRECISION/OCG/lizihao/projects/R/library")
library(JASPAR2020)
### jaspar function
getJasparMotifs2020 <- function(species = "Homo sapiens", 
                              collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}
###

.libPaths("/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/anaconda3/lib/R/library")
data <- readRDS(paste0("/hwfssz5/ST_PRECISION/OCG/lizihao/projects/brain/ATAC/PeakMatrix/",sub,"PeakMatrix"))
data <- data@assays@data$PeakMatrix
peakset <- readRDS("/hwfssz5/ST_PRECISION/OCG/lizihao/projects/brain/ATAC/PeakMatrix/peakset_character.RDS")
rownames(data) <- peakset
#rownames(data) <- gsub("chr","MFA",rownames(data))
peakfile <- "/hwfssz5/ST_PRECISION/OCG/lizihao/projects/brain/ATAC/PeakMatrix/peak.bed"
peaks <- getPeaks(peakfile, sort_peaks = TRUE)
#data <- as.matrix(data)
fragment_counts <- SummarizedExperiment(assays = 
                                          list(counts = data),
                                        rowRanges = peaks)


library(BSgenome.Mfascicularis.NCBI.5.0)
example_counts <- addGCBias(fragment_counts, 
                            genome = BSgenome.Mfascicularis.NCBI.5.0)
counts_filtered <- filterPeaks(example_counts, non_overlapping = TRUE)

motifs <- getJasparMotifs2020()
library(motifmatchr)
motif_ix <- matchMotifs(motifs, counts_filtered, 
                        genome = BSgenome.Mfascicularis.NCBI.5.0)
dev <- computeDeviations(object = counts_filtered, 
                                 annotations = motif_ix)
variability <- computeVariability(dev)
saveRDS(variability,paste0("/hwfssz5/ST_PRECISION/OCG/lizihao/projects/brain/ATAC/PeakMatrix/chromvar2020/",sub,"var.RDS"))
saveRDS(deviationScores(dev),paste0("/hwfssz5/ST_PRECISION/OCG/lizihao/projects/brain/ATAC/PeakMatrix/chromvar2020/",sub,"deviationScores.RDS"))
