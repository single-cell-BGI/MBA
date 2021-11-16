library(ArchR)
projHeme2 <- readRDS("Save-ArchR-Project.rds")
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI_15000", 
    iterations = 10, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 15000, 
    dimsToUse = 1:30
)
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI_15000",
    name = "Harmony_15000",
    groupBy = c("Monkey","Region")
)
projHeme2 <- addClusters(
     input = projHeme2,
     reducedDims = "Harmony_15000",
     method = "Seurat",
     name = "Clusters_15000",
     resolution = 0.8
)
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony_15000",
    name = "UMAPHarmony_15000",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)

