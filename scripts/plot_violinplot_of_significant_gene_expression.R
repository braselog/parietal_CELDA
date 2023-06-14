# plot significant genes
library(celda)
library(ggplot2)
library(cowplot)
library(qs)
library(Seurat)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
# test
# args <- c("deleteme.txt", "astro", "~/SingleCellProjects/MyProjects/SexDifferences/celdaAnalyses/Astro/celda_astro_splitCell.qs", 65, 8)

hits = NULL
try({
    hits <- read.table(args[1], sep = ",", header = T)
})
if (is.null(hits)) {
    write("", args[6])
} else {
    gene_clust <- colnames(hits[, as.numeric(hits["BHcombined", ]) < 0.05])
    m <- readRDS(sprintf("~/SingleCellProjects/dataObjects/%s.rds", args[2]))
    sce <- qread(args[3], nthreads = 8)
    L <- read.table(file = args[4], col.names = "L")$L
    K <- read.table(file = args[5], col.names = "K")$K
    k <- subsetCeldaList(sce, list(L = L, K = K))
    rm(sce)
    m$celdaClusters <- paste0("K", celdaClusters(k))
    rm(k)
    m <- subset(m, subset = Sample_ID %in% names(table(m@meta.data$Sample_ID)[table(m@meta.data$Sample_ID) > 60])) # 60]))
    m$Clusters <- Idents(m)

    df <- data.frame(t(GetAssayData(m, slot = "counts")[unique(unlist(hits["Gene", gene_clust])), ]), CeldaCluster = m$celdaClusters, Sex = as.character(m$Gender), Status = unlist(lapply(m$Status, function(x) strsplit(as.character(x), "_")[[1]][2])), check.names = F)
    df$Status[df$Status == "AD"] <- "sAD"
    df$Sex[df$Sex == 1] <- "Male"
    df$Sex[df$Sex == 2] <- "Female"
    df$Status_Sex <- paste(df$Status, df$Sex, sep = "-")
    df$Status_Sex <- factor(df$Status_Sex, levels = c("CO-Female", "sAD-Female", "CO-Male", "sAD-Male"))
    # colnames(df)[which(colnames(df) == "FAM189A2")] <- gene_clust[which(gene_clust == "FAM189A2")] <- "ENTREP1"
    plts <- lapply(gene_clust, function(gc) {
        df_tmp <- df[df$CeldaCluster == hits["CeldaCluster", gc], ] # filter the data to only include the cluster of interest
        p <- ggplot(na.omit(df_tmp), aes(x = Status_Sex, y = log(UQ(as.name(hits["Gene", gc])) + 1), fill = Sex)) +
            geom_jitter(size = 0.5) +
            geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
            stat_summary(fun = "mean", geom = "crossbar", width = 0.75, color = "red") +
            ggtitle(gc) +
            theme_minimal() +
            xlab("Status-Sex") +
            theme(axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        return(p)
    })
    q <- plot_grid(plotlist = plts, ncol = 3, nrow = ceiling(length(gene_clust) / 3))
    ggsave(q, file = args[6], w = 12, h = 4.5 * ceiling(length(gene_clust) / 3), limitsize = F)
}
