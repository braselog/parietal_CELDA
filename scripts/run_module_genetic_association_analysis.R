library(qs)
library(lme4)
library(lmerTest)
library(celda)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
# test
# args <- c("results/astro/K_celda_astro.qs", "results/astro/L_astro.txt", "results/astro/K_astro.txt", "results/astro/celda_gene_module_status_model_celdaCluster_astro.qs", "results/astro/celda_gene_module_status_model_celdaCluster_astro.csv")

sce <- qread(args[1], nthreads = 8)
L <- read.table(file = args[2], col.names = "L")$L
K <- read.table(file = args[3], col.names = "K")$K
k <- subsetCeldaList(sce, list(L = L, K = K))
rm(sce)

mdata <- readRDS("~/SingleCellProjects/dataObjects/allCellTypes_garnettCleaned_finalObj_294114nuclei_metaData.rds")
mdata <- mdata[colnames(k), ]
fm <- factorizeMatrix(k)
allModDF <- data.frame(t(fm$proportions$cell), Sex = mdata$Gender, celdaCluster = paste0("K", celdaClusters(k)), cellType = mdata$cellType, cellState = mdata$cellState, Status = mdata$Status, Sample = mdata$Sample_ID, rs1582763 = mdata$MS4, APOE4 = as.numeric(grepl("4", mdata$nAPOE)))
allModFilt <- allModDF
allModFilt$Status <- factor(allModFilt$Status, levels = c("Neuro_CO", "Neuro_AD", "Neuro_ADAD", "Neuro_Presympt", "Neuro_OT"))
allModFilt$rs1582763 <- as.numeric(factor(allModFilt$rs1582763, levels = c("GG", "AG", "AA")))

# FIXME: delete this chunk. we like the celdaCluster model better
## the following code is to run the sex*AD model for each module in all astrocytes
# myHits <- list()
# tmp <- lapply(seq_len(args[2]), function(idx) {
#    mod <- paste0("L", idx)
#    re_ln2 <- lmer(paste0(mod, " ~ Sex + Status + rs1582763 + APOE4 + (1|Sample)"), data = allModFilt, REML = TRUE)
#    z <- summary(re_ln2)
#    myHits[[mod]] <<- z$coefficients[-1, ][z$coefficients[-1, "Pr(>|t|)"] < 0.05, , drop = F]
#    return(z$coefficients)
# })
# qsave(tmp, file = "~/SingleCellProjects/MyProjects/SexDifferences/celdaAnalyses/Astro/celda_gene_module_status_model.qs", nthreads = 8)

# the following code is to run the sex*AD model for each module in all astrocytes accounting for celda cluster
myHits <- list()
i <- 0
tmp <- lapply(seq_len(L), function(idx) {
    i <<- i + 1
    mod <- paste0("L", idx)
    re_ln2 <- lmer(paste0(mod, " ~ Sex + celdaCluster + Status + rs1582763 + APOE4 + (1|Sample)"), data = allModFilt, REML = TRUE)
    z <- summary(re_ln2)$coefficients
    z2 <- data.frame(z[-grep("celdaClusterK", rownames(z)), ])
    z2$mod <- mod
    z2 <- cbind(rownames(z2), z2)
    myHits[[mod]] <<- z2[-1, ][z2[-1, "Pr...t.."] < 0.05, , drop = F]
    return(z)
})
myHits2 <- do.call(rbind, myHits)
qsave(tmp, file = args[4], nthreads = 8)
write.table(myHits2, file = args[5], sep = ",", quote = F)


# TODO: do i need this chunk? i have not included it in the snakemake file yet.
## plot the L17 module score by Status using ggplot2
# library(ggplot2)
# p <- ggplot(allModFilt, aes(x = Status, y = L17)) +
#    geom_jitter(width = 0.2, height = 0.01, alpha = 0.5) +
#    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#    ggtitle("L17 module score by Status")
## save plot
# ggsave(p, file = "~/SingleCellProjects/MyProjects/SexDifferences/celdaAnalyses/Astro/celdaModule_L17_by_status.pdf", w = 7, h = 7)
