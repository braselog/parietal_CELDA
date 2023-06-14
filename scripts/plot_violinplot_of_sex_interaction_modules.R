library(celda)
library(ggplot2)
library(cowplot)
library(qs)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
# test
# args <- c("results/astro/K_celda_astro.qs", "results/astro/L_astro.txt", "results/astro/K_astro.txt", "results/astro/sig_module_cluster_sAD_pvals_astro.qs", "results/astro/violinplot_of_sex_interaction_modules_astro.pdf")
# setwd("celdaAnalyses/DiscoveryPipeline/")

sce <- qread(args[1], nthreads = 8)
L <- read.table(file = args[2], col.names = "L")$L
K <- read.table(file = args[3], col.names = "K")$K
k <- subsetCeldaList(sce, list(L = L, K = K))
rm(sce)
fm <- factorizeMatrix(k)

# read in the sAD interaction pvalues for the modules
sAD_pvals <- qread(args[4])
if (sum(sAD_pvals < 0.05) == 0) {
    write("", args[5])
} else {
    mods <- unlist(lapply(which(sAD_pvals < 0.05), function(idx) {
        return(strsplit(names(sAD_pvals)[idx], "_")[[1]][1])
    }))
    cStates <- unlist(lapply(which(sAD_pvals < 0.05), function(idx) {
        return(strsplit(names(sAD_pvals)[idx], "_")[[1]][2])
    }))
    # these are interesting to look at
    # apply(fm$posterior$module, 2, function(mod) {
    #    return(sum(mod >= 0.01))
    # })[mods] # how many genes contribute more than 1% to the total module score
    # apply(fm$posterior$module, 2, function(mod) {
    #    return(sum(mod[mod >= 0.01]))
    # })[mods] # using those genes what is their combined contribution
    mdata <- readRDS("~/SingleCellProjects/dataObjects/allCellTypes_garnettCleaned_finalObj_294114nuclei_metaData.rds")
    mdata <- mdata[colnames(k), ]
    allModFilt <- data.frame(t(fm$proportions$cell), Sex = mdata$Gender, celdaCluster = paste0("K", celdaClusters(k)), cellType = mdata$cellType, cellState = mdata$cellState, Status = mdata$Status, Sample = mdata$Sample_ID)

    # violin plots of these hits
    plts <- lapply(1:length(mods), function(ind) {
        df <- data.frame(Module_Score = allModFilt[, mods[ind]], Status_Sex = paste(allModFilt$Status, allModFilt$Sex, sep = "_"), Sex = allModFilt$Sex, Status = allModFilt$Status, celdaCluster = allModFilt$celdaCluster, cellState = allModFilt$cellState, Sample = allModFilt$Sample)
        df <- df[df$celdaCluster == cStates[ind], ]
        df$Status_Sex <- factor(df$Status_Sex, levels = c("Neuro_CO_1", "Neuro_AD_1", "Neuro_ADAD_1", "Neuro_CO_2", "Neuro_AD_2", "Neuro_ADAD_2"))
        p <- ggplot(na.omit(df), aes(x = Status_Sex, y = Module_Score, fill = Sex)) +
            geom_jitter() +
            geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.75) +
            stat_summary(fun = "mean", geom = "crossbar", width = 0.75, color = "red") +
            ggtitle(sprintf("%s : %s", mods[ind], cStates[ind]))
        return(p)
    })
    q <- plot_grid(plotlist = plts, ncol = 3, nrow = ceiling(length(mods) / 3))
    ggsave(q, file = args[5], w = 21, h = 7 * ceiling(length(mods) / 3), limitsize = FALSE)
}
