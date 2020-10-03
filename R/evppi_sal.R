## Code taken from BCEA package
## Baio, G., Berardi, A., & Heath, A. (2017). Bayesian cost-effectiveness analysis with the R package BCEA. New York: Springer.
## https://github.com/giabaio/BCEA

evppi_sal <- function(outputs, inputs, pars, ...){
    n.seps <- list(...)$n.seps
    if (is.null(n.seps)) n.seps <- 1

    U <- form_nbarray(outputs, inputs)
    nsim <- dim(U)[1]
    nk <- dim(U)[2]
    nopt <- dim(U)[3]

    if (length(pars) > 1)
        stop("`method=\"sal\" only works for single-parameter EVPPI")
    
    param <- inputs[, pars]
    o <- order(param)
    param <- param[o]
    nSegs <- matrix(1, nopt, nopt)
    nSegs[1, 2] <- n.seps
    nSegs[2, 1] <- n.seps
    res <- segPoints <- numeric()
    for (k in 1:nk) {
        nbs <- U[, k, ]
        nbs <- nbs[o, ]
        for (i in 1:(nopt - 1)) {
            for (j in (i + 1):nopt) {
                cm <- cumsum(nbs[, i] - nbs[, j])/nsim
                if (nSegs[i, j] == 1) {
                    l <- which.min(cm)
                    u <- which.max(cm)
                    if (cm[u] - max(cm[1], cm[nsim]) > min(cm[1],
                                                           cm[nsim]) - cm[l]) {
                        segPoint <- u
                    }
                    else {
                        segPoint <- l
                    }
                    if (segPoint > 1 && segPoint < nsim) {
                        segPoints <- c(segPoints, segPoint)
                    }
                }
                if (nSegs[i, j] == 2) {
                    distMaxMin <- 0
                    distMinMax <- 0
                    minL <- Inf
                    maxL <- -Inf
                    for (sims in 1:nsim) {
                        if (cm[sims] > maxL) {
                            maxLP <- sims
                            maxL <- cm[sims]
                        }
                        else {
                            if (maxL - cm[sims] > distMaxMin) {
                                distMaxMin <- maxL - cm[sims]
                                segMaxMinL <- maxLP
                                segMaxMinR <- sims
                            }
                        }
                        if (cm[sims] < minL) {
                            minLP <- sims
                            minL <- cm[sims]
                        }
                        else {
                            if (cm[sims] - minL > distMinMax) {
                                distMinMax <- cm[sims] - minL
                                segMinMaxL <- minLP
                                segMinMaxR <- sims
                            }
                        }
                    }
                    siMaxMin <- cm[segMaxMinL] + distMaxMin +
                        (cm[nsim] - cm[segMaxMinR])
                    siMinMax <- -cm[segMaxMinL] + distMinMax -
                        (cm[nsim] - cm[segMinMaxR])
                    if (siMaxMin > siMinMax) {
                        segPoint <- c(segMaxMinL, segMaxMinR)
                    }
                    else {
                        segPoint <- c(segMinMaxL, segMinMaxR)
                    }
                    if (segPoint[1] > 1 && segPoint[1] <
                        nsim) {
                        segPoints <- c(segPoints, segPoint[1])
                    }
                    if (segPoint[2] > 1 && segPoint[2] <
                        nsim) {
                        segPoints <- c(segPoints, segPoint[2])
                    }
                }
            }
        }
        if (length(segPoints) > 0) {
            segPoints2 <- unique(c(0, segPoints[order(segPoints)],
                                   nsim))
            res[k] <- 0
            for (j in 1:(length(segPoints2) - 1)) {
                res[k] <- res[k] + max(colSums(matrix(nbs[(1 +
                                                           segPoints2[j]):segPoints2[j + 1], ],
                                                      ncol = nopt)))/nsim
            }
            res[k] <- res[k] - max(colMeans(nbs))
        }
        else {
            res[k] <- 0
        }
    }
    res
}
