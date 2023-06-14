library(qs)
library(lme4)
library(lmerTest)
library(celda)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
# test
# setwd("celdaAnalyses/DiscoveryPipeline/")
# args <- c("results/astro/K_celda_astro.qs", "results/astro/L_astro.txt", "results/astro/K_astro.txt", "results/astro/all_module_cluster_results_astro.txt", "results/astro/sig_module_cluster_results_astro.txt", "results/astro/sig_module_cluster_sAD_pvals_astro.qs")

sce <- qread(args[1], nthreads = 8)
L <- read.table(file = args[2], col.names = "L")$L
K <- read.table(file = args[3], col.names = "K")$K
k <- subsetCeldaList(sce, list(L = L, K = K))
rm(sce)

mdata <- readRDS("~/SingleCellProjects/dataObjects/allCellTypes_garnettCleaned_finalObj_294114nuclei_metaData.rds")
mdata <- mdata[colnames(k), ]
fm <- factorizeMatrix(k)
allModDF <- data.frame(t(fm$proportions$cell), Sex = mdata$Gender, celdaCluster = paste0("K", celdaClusters(k)), cellType = mdata$cellType, cellState = mdata$cellState, Status = mdata$Status, Sample = mdata$Sample_ID)
allModFilt <- allModDF
allModFilt$Status <- factor(allModFilt$Status, levels = c("Neuro_CO", "Neuro_AD", "Neuro_ADAD", "Neuro_Presympt", "Neuro_OT"))

# determine the module-cluster pairs to test
pairs <- apply(fm$proportions$cellPopulation, 1, function(r) {
    z <- kmeans(r, centers = 2)
    return(names(z$cluster)[which(z$cluster == names(table(z$cluster))[which(table(z$cluster) == min(table(z$cluster)))])]) # get the cluster names of the clusters which have high expression of the module
})
clusts <- pairs

# the following code is to run the sex*AD model for each module
allModuleResults <- lapply(seq_len(nrow(fm$proportions$cellPopulation)), function(idx) {
    mod <- paste0("L", idx)
    cs <- clusts[[idx]]
    results <- lapply(cs, function(cstate) {
        re_ln2 <- NA
        try({
            re_ln2 <- lmer(paste0(mod, " ~ Sex + Status +  Sex*Status + (1|Sample)"), data = allModFilt[allModFilt$celdaCluster == cstate, ], REML = TRUE)
        }) # nebula(m$DAM_up1,m@meta.data$Sample_ID,pred=pred,method='LN',model='gaussain')#'NBLMM')
        return(summary(re_ln2))
    })
    names(results) <- cs
    return(results)
})
names(allModuleResults) <- paste0("L", seq_len(nrow(fm$proportions$cellPopulation)))
qsave(allModuleResults, file = args[4], nthreads = 8)

# the following code prints the results that are less than 0.05
sAD_pvals <- c()
write("", file = args[5])
tmp <- lapply(names(allModuleResults), function(mod) {
    lapply(names(allModuleResults[[mod]]), function(res) {
        c <- NA
        try({
            c <- allModuleResults[[mod]][[res]]$coefficients
        })
        if (!sum(is.na(c))) {
            if ("Sex2:StatusNeuro_AD" %in% rownames(c)) {
                sAD_pvals <<- c(sAD_pvals, c["Sex2:StatusNeuro_AD", "Pr(>|t|)"])
                names(sAD_pvals)[length(sAD_pvals)] <<- sprintf("%s_%s", mod, res)
                if (sum(c["Sex2:StatusNeuro_AD", "Pr(>|t|)"] < 0.05) > 0) {
                    write(sprintf("%s_%s", mod, res), file = args[5], append = T)
                    write.table(c[c("Sex2", "Sex2:StatusNeuro_AD"), ], file = args[5], append = T, quote = F, sep = "\t", row.names = T, col.names = T)
                }
            }
        }
    })
})
qsave(sAD_pvals, file = args[6], nthreads = 8)
