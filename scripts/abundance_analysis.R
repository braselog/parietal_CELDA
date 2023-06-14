# use propeller https://github.com/Oshlack/speckle

library(speckle)
library(qs)
library(celda)
library(limma)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

sce <- qread(args[1], nthreads = 8)
L <- read.table(file = args[2], col.names = "L")$L
K <- read.table(file = args[3], col.names = "K")$K
k <- subsetCeldaList(sce, list(L = L, K = K))
rm(sce)

k$cluster <- celdaClusters(k)
k$sample <- k$Sample_ID
k$group <- k$Gender

props <- getTransformedProps(k$cluster, k$sample, transform = "asin")
props$TransformedProps <- props$Proportions^(1 / 3)
groups <- unique(cbind(k$group, k$sample))
rownames(groups) <- groups[, 2]
groups <- groups[colnames(props$TransformedProps), 1]
groups <- factor(groups, levels = c(1, 2))
Status <- unique(cbind(k$Status, k$sample))
rownames(Status) <- Status[, 2]
Status <- Status[colnames(props$TransformedProps), 1]
Status <- factor(Status, levels = c("Neuro_CO", "Neuro_Presympt", "Neuro_AD", "Neuro_ADAD", "Neuro_OT"))
design <- model.matrix(~ 0 + groups + Status)
mycontr <- makeContrasts(groups1 - groups2, levels = design)

results <- propeller.ttest(prop.list = props, design = design, contrasts = mycontr, robust = T, trend = T, sort = T)
write.table(results, file = args[4], sep = ",", quote = F, row.names = T, col.names = T)
