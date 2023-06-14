library(qs)
library(celda)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
# test
# args <- c("~/SingleCellProjects/MyProjects/SexDifferences/celdaAnalyses/Astro/celda_astro_splitCell.qs", 65, 8, "deleteme.txt", "astro")

sce <- qread(args[1], nthreads = 8)
L <- read.table(file = args[2], col.names = "L")$L
K <- read.table(file = args[3], col.names = "K")$K
k <- subsetCeldaList(sce, list(L = L, K = K))
rm(sce)
fm <- factorizeMatrix(k)

# read in the significant module cluster pairs
MOD_CLUST <- unlist(lapply(readLines(args[4]), function(l) if (startsWith(l, "L")) l else NULL))
if (is.null(MOD_CLUST)) {
    write("", file = args[6])
} else {
    MOD_HITS <- unlist(lapply(MOD_CLUST, function(l) strsplit(l, "_")[[1]][1]))
    CLUSTER_HITS <- unlist(lapply(MOD_CLUST, function(l) strsplit(l, "_")[[1]][2]))
    names(MOD_HITS) <- MOD_CLUST
    names(CLUSTER_HITS) <- MOD_CLUST


    modListFull <- lapply(colnames(fm$posterior$module), function(mod) {
        return(fm$posterior$module[, mod])
    })
    modListFull <- lapply(modListFull, function(mod) {
        return(mod[order(mod, decreasing = T)])
    })
    names(modListFull) <- colnames(fm$posterior$module)
    mod80Full <- lapply(modListFull, function(mod) { # using this filtering method: on average 23% of the gene for each module make up 80% of the overall signal
        ind <- 1
        while (sum(mod[1:ind]) < 0.8) ind <- ind + 1
        return(mod[1:ind])
    })
    names(mod80Full) <- colnames(fm$posterior$module)
    genes <- lapply(MOD_HITS, function(mod) {
        return(names(mod80Full[[mod]]))
    })
    names(genes) <- MOD_CLUST
    rm(modListFull)
    rm(mod80Full)

    library(nebula)
    library(Seurat)
    library(qs)

    pvalue.extreme <- function(z) {
        log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
        log10.pvalue <- log.pvalue / log(10) ## from natural log to log10
        mantissa <- 10^(log10.pvalue %% 1)
        # mantissa[!is.finite(mantissa)] <- 1
        exponent <- log10.pvalue %/% 1
        # exponent[!is.finite(exponent)] <- 0
        ## or return(c(mantissa,exponent))
        # return(sprintf("p value is %1.2f times 10^(%d)",mantissa,exponent))
        return(sprintf("%1.2fE%d", mantissa, exponent))
    }

    BH.extreme <- function(Zs) {
        len <- length(Zs)
        BHs <- unlist(lapply(Zs, function(z) {
            log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) + log(len) - log(which(Zs == z))
            log10.pvalue <- log.pvalue / log(10) ## from natural log to log10
            mantissa <- 10^(log10.pvalue %% 1)
            exponent <- log10.pvalue %/% 1
            ## or return(c(mantissa,exponent))
            # return(sprintf("p value is %1.2f times 10^(%d)",mantissa,exponent))
            return(sprintf("%1.2fE%d", mantissa, exponent))
        }))
        for (i in 1:length(BHs)) {
            change <- which(as.numeric(BHs[i]) < as.numeric(BHs[1:(i - 1)]))
            if (length(change) > 0) {
                BHs[change] <- BHs[i]
            }
        }
        return(BHs)
    }

    m <- readRDS(sprintf("~/SingleCellProjects/dataObjects/%s.rds", args[5]))
    m$celdaClusters <- paste0("K", celdaClusters(k))
    rm(k)
    m <- subset(m, subset = Sample_ID %in% names(table(m@meta.data$Sample_ID)[table(m@meta.data$Sample_ID) > 60])) # 60]))
    m$Clusters <- Idents(m)

    genes <- lapply(genes, function(g) {
        return(g[g %in% rownames(m)])
    })

    m$Sex <- as.factor(m$Gender)
    m$AOD <- scale(as.numeric(m$AOD))
    m$PMI <- scale(as.numeric(m$PMI))
    m$Status <- factor(m$Status, levels = c("Neuro_CO", "Neuro_AD", "Neuro_ADAD", "Neuro_Presympt", "Neuro_OT"))

    # run the linear model for each gene within the associated cluster
    allDE <- NULL
    msub <- NULL
    tmp <- lapply(MOD_CLUST, function(L_K) {
        msub <<- subset(m, subset = celdaClusters == CLUSTER_HITS[L_K])
        statCount <- table(msub$Status)
        if (sum(statCount < 50) > 0) msub <<- subset(msub, subset = Status %in% names(statCount)[statCount >= 50])
        pred <- model.matrix(~ Sex + Status + Sex * Status, data = msub@meta.data)
        pred <- pred[, colSums(pred) != 0] # remove predictors with 0 variance
        re_ln2 <- nebula(GetAssayData(msub, slot = "counts")[genes[[L_K]], ], msub@meta.data$Sample_ID, pred = pred, method = "LN", model = "NBLMM")
        tmp2 <- data.frame("logFC_Sex.StatusNeuro_AD" = re_ln2$summary[, "logFC_Sex2:StatusNeuro_AD"])
        tmp2$z.score <- re_ln2$summary[, "logFC_Sex2:StatusNeuro_AD"] / re_ln2$summary[, "se_Sex2:StatusNeuro_AD"]
        tmp2$z.score[!is.finite(tmp2$z.score)] <- 0.0
        tmp2 <- tmp2[order(abs(tmp2$z.score), decreasing = T), ]
        tmp2$p.value <- pvalue.extreme(tmp2$z.score) # re_ln2$summary[, "p_Sex2:StatusNeuro_AD"] #
        tmp2$BH <- BH.extreme(tmp2$z.score) # p.adjust(re_ln2$summary[, "p_Sex2:StatusNeuro_AD"], "BH")
        merged2 <- merge(tmp2, re_ln2$summary, by.x = "logFC_Sex.StatusNeuro_AD", by.y = "logFC_Sex2:StatusNeuro_AD")
        rownames(merged2) <- paste(merged2$gene, CLUSTER_HITS[L_K], sep = "_")
        merged2 <- merged2[order(abs(merged2$z.score), decreasing = T), ]
        merged2$CeldaCluster <- CLUSTER_HITS[L_K]
        merged2$GeneModule <- MOD_HITS[L_K]
        colnames(merged2)[colnames(merged2) == "gene"] <- "Gene"
        allDE <<- rbind(allDE, merged2[, c("Gene", "CeldaCluster", "GeneModule", "logFC_Sex2", "logFC_StatusNeuro_AD", "logFC_Sex.StatusNeuro_AD", "p_Sex2", "p_StatusNeuro_AD", "p.value", "BH")])
        return(merged2)
    })
    allDE$BHcombined <- p.adjust(allDE$p.value, "BH")
    allDE <- allDE[order(allDE$BHcombined, decreasing = F), ]
    write.table(t(allDE[seq(sum(as.numeric(allDE$p.value) < 0.05)), ]), file = args[6], sep = ",", quote = F, row.names = T, col.names = T)
}
