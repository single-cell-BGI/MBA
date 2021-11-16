library(ArchR)
library(Seurat)

proj <- loadArchRProject(subATAC))
rna <- readRDS(subRNA)
projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI_15000",
    seRNA = rna,
    addToArrow = TRUE,
    groupRNA = "celltype_5",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    threads = 1,
    force = TRUE,
)
projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI_15000"
)
p2g <- getPeak2GeneLinks(
    ArchRProj = projHeme5,
    corCutOff = 0.2,
    resolution = 1,
    returnLoops = FALSE
)
p2geneDF <- p2g
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
sub_peak <- readRDS("/PeakCalls/sub_peak_path-reproduciblePeaks.gr.rds")
sub_peak_1k_add <- sub_peak[sub_peak$distToTSS >= 1000,]
sub_peak_1k_low <- sub_peak[sub_peak$distToTSS < 1000,]
sub_peak_1k_add <- as.data.frame(sub_peak_1k_add)
sub_peak_1k_low <- as.data.frame(sub_peak_1k_low)
sub_peak <- as.data.frame(sub_peak)
sub_peak_1k_add$peakName <- paste0(sub_peak_1k_add$seqnames,"_",sub_peak_1k_add$start,"_",sub_peak_1k_add$end)
sub_peak_1k_low$peakName <- paste0(sub_peak_1k_low$seqnames,"_",sub_peak_1k_low$start,"_",sub_peak_1k_low$end)

sub_peak$peakName <- paste0(sub_peak$seqnames,"_",sub_peak$start,"_",sub_peak$end)
p2geneDF_sub <- p2geneDF[p2geneDF$peakName %in% sub_peak$peakName,]
p2geneDF_sub[p2geneDF_sub$peakName %in% sub_peak_1k_add$peakName,"is_promoter"] <- "F"
p2geneDF_sub[p2geneDF_sub$peakName %in% sub_peak_1k_low$peakName,"is_promoter"] <- "T"

