library(celda)
library(qs)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

moduleSplit <- qread(args[1], nthreads = 8)
L <- read.table(file = args[2], col.names = "L")$L
sce <- subsetCeldaList(moduleSplit, list(L = L))
cModules <- celdaModules(sce)
rm(moduleSplit)
sce <- recursiveSplitCell(sce, initialK = 3, maxK = 25, sampleLabel = sce$Sample_ID, yInit = cModules)
gc()
qsave(sce, file = args[3], nthreads = 8)
