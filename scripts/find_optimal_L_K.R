library(qs)
library(celda)

## elbow functions from https://rdrr.io/bioc/celda/src/R/elbow.R
.dist2d <- function(a, b, c) {
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1, v2)
    d <- abs(det(m)) / sqrt(sum(v1 * v1))
    return(d)
}

.secondDerivativeEstimate <- function(v) {
    nv <- length(v)
    res <- rep(NA, nv)
    for (i in seq(2, nv - 1)) {
        res[i] <- v[i + 1] + v[i - 1] - (2 * v[i])
    }
    return(res)
}

.curveElbow <- function(var, perplexity, pvalCutoff = 0.05) {
    len <- length(perplexity)
    a <- c(var[1], perplexity[1])
    b <- c(var[len], perplexity[len])
    res <- rep(NA, len)
    for (i in seq_along(var)) {
        res[i] <- .dist2d(c(var[i], perplexity[i]), a, b)
    }
    elbow <- which.max(res)
    ix <- var > var[elbow]
    perplexitySde <- .secondDerivativeEstimate(perplexity)
    perplexitySdeSd <- stats::sd(perplexitySde[ix], na.rm = TRUE)
    perplexitySdeMean <- mean(perplexitySde[ix], na.rm = TRUE)
    perplexitySdePval <-
        stats::pnorm(perplexitySde,
            mean = perplexitySdeMean,
            sd = perplexitySdeSd,
            lower.tail = FALSE
        )
    # other <- which(ix & perplexitySdePval < pvalCutoff)
    return(elbow = var[elbow])
}

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

sce <- qread(args[1], nthreads = 8)

pdf(args[2], width = 15)
plotGridSearchPerplexity(sce)
plotRPC(sce)
dev.off()

# find the inflection point in the rate of change of perplexity
p <- plotGridSearchPerplexity(sce)
pc <- diff(p$data$perplexity) # rate of change of perplexity
x <- seq_len(length(pc))
y <- pc
# plot(x, y, type = "l")
lo <- loess(y ~ x)
# xl <- seq(min(x), max(x), (max(x) - min(x)) / 1000)
out <- predict(lo, x)
# lines(xl, out, col = "red", lwd = 2)

# Find inflection point
infl <- c(FALSE, diff(diff(out) > 0) != 0)
# only grab the first inflection point (the one closest to the elbow)
infl <- as.numeric(as.character(p$data$L[round(xl[which(infl)])]))[1]
if (is.na(infl) && !grepl("K", args[1])) {
    cat("ERROR: Inflection point not found!!!\n")
    # infl <- max(as.numeric(as.character(p$data$L)))
}

# find the elbow point in the perplexity
elbow <- round(.curveElbow(x, out))
# elbow <- round(.curveElbow(x, pc))
# elbow <- round(.curveElbow(as.numeric(as.character(p$data$K)), p$data$perplexity))
if (grepl("K", args[1])) {
    write(p$data$K[elbow], file = args[3]) # for the number of clusters use the elbow to prevent overclustering
} else {
    L <- round((infl + as.numeric(as.character(p$data$L[elbow]))) / 2)
    write(L, file = args[3]) # for the number of modules use the average of the elbow and inflection point
}
