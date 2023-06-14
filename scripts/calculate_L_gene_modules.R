# Load libraries
library(Seurat)
library(celda)
library(qs)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

# FIXME: I need to filter the samples by those that are european ancestry
sce <- readRDS(paste0("~/SingleCellProjects/dataObjects/", args[1], ".rds"))
sce <- as.SingleCellExperiment(sce, assay = "SCT")
#
useAssay <- "counts"
altExpName <- "featureSubset"
sce <- selectFeatures(sce, minCount = 3, minCell = 3, useAssay = useAssay, altExpName = altExpName)
nrow(altExp(sce, altExpName))

moduleSplit <- recursiveSplitModule(sce, useAssay = useAssay, altExpName = altExpName, initialL = 10, maxL = 125, sampleLabel = sce$Sample_ID, verbose = T, zInit = as.numeric(factor(sce$seurat_clusters)))
qsave(moduleSplit, file = args[2], nthreads = 8)
