library(ArchR)
args=commandArgs(T)
sub <- args[1]
proj4 = loadArchRProject(paste0("/hwfssz1/ST_PRECISION_COVID/3.lizihao/brain/ATAC/",sub,"Save-ProjMonkey-celltype-match"))
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj4, 
    useMatrix = "PeakMatrix", 
    groupBy = "rcelltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks,paste0("/hwfssz1/ST_PRECISION_COVID/3.lizihao/brain/ATAC/",sub,"Save-ProjMonkey-celltype-match/",sub,"markersPeaks.RDS"))
markerList <- getMarkers(markersPeaks, cutOff = "Pval <= 0.01 & Log2FC >= 1")
saveRDS(markerList,paste0("/hwfssz1/ST_PRECISION_COVID/3.lizihao/brain/ATAC/",sub,"Save-ProjMonkey-celltype-match/",sub,"markerList_p0.01fc1.RDS"))
