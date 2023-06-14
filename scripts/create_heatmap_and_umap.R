library(celda)
library(qs)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

sce <- qread(args[1], nthreads = 8)
L <- read.table(file = args[2], col.names = "L")$L
K <- read.table(file = args[3], col.names = "K")$K
k <- subsetCeldaList(sce, list(L = L, K = K))
rm(sce)

k <- celdaUmap(k)

pdf(args[4])
plotDimReduceCluster(k, reducedDimName = "celda_UMAP", labelClusters = TRUE)
dev.off()

pdf(args[5], height = 12)
celdaProbabilityMap(k)
dev.off()

write.table(table(k$seurat_clusters, celdaClusters(k)), file = args[6], sep = ",", quote = FALSE)
write("", file = args[6], append = TRUE, sep = ",")
write.table(table(k$Sample_ID, celdaClusters(k)), file = args[6], append = T, sep = ",", quote = FALSE)
