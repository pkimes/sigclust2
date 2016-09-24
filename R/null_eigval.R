#' Eigenvalue estimation for null Gaussian based testing procedures
#'
#' Function to compute the eigenvalues of the null Gaussian distribution for
#' significance of clustering testing procedures which rely on a null Gaussian
#' factor model assumption. When the number of observations is substantially
#' greater than the number of features, the sample covariance matrix should be
#' used.
#' 
#' @param x a matrix of size n by p containing the original data.
#' @param n an integer number of samples.
#' @param p an integer number of features/covariates.
#' @param icovest an integer between 1 and 3 corresponding to the covariance
#'        estimation procedure to use. See details for more
#'        information on the possible estimation procedures. (default = 1)
#' @param bkgd_pca a logical value specifying whether to use scaled PCA scores
#'        over raw data to estimate the background noise. When FALSE, raw estimate
#'        is used; when TRUE, minimum of PCA and raw estimates is used.
#'        (default = FALSE)
#' 
#' @return
#' The function returns a list of estimated parameters for the null Gaussian
#' distribution used in significance of clustering testing. The list includes:
#' \itemize{
#'     \item \code{eigval_dat}: eigenvalues for sample covariance matrix
#'     \item \code{backvar}: background noise, sigma_b^2
#'     \item \code{eigval_sim}: eigenvalues to be used for simulation
#' }
#' 
#' @details
#' The following possible options are given for null covariance estimation
#' \enumerate{
#'     \item soft thresholding: recommended approach described in Huang et al. 2014
#'     \item sample: uses sample covariance matrix, equivalent to soft and hard options
#'           when n > p, but when p > n, will produce conservative results, i.e. less
#'           significant p-values
#'     \item hard thresholding: approach described in Liu et al. 2008, no longer
#'           recommended - retained for historical purposes
#' }
#' 
#' @references
#' \itemize{
#'     \item Huang, H., Liu, Y., Yuan, M., and Marron, J. S. (2014).
#'           Statistical Significance of Clustering using Soft Thresholding.
#'           Journal of Computational and Graphical Statistics, preprint.
#'     \item Liu, Y., Hayes, D. N., Nobel, A. B., and Marron, J. S. (2008).
#'           Statistical Significance of Clustering for High-Dimension, Low-Sample Size Data.
#'           Journal of the American Statistical Association, 103(483):1281-1293.
#' }
#'
#' @export
#' @author Patrick Kimes
null_eigval <- function(x, n, p, icovest = 1, bkgd_pca = FALSE) {

    if (!(icovest %in% 1:3)) {
        warning("icovest should be 1, 2 or 3. Using default value: 1.")
        icovest <- 1
    }
    
    if (nrow(x) != n | ncol(x) != p)
        stop("Wrong size of matrix x!")
    
    ## compute background based on raw data
    ## or min of raw data and pca scores
    mad1 <- mad(as.matrix(x))
    if (bkgd_pca) {
        mad1 <- min(mad1, mad(as.matrix(prcomp(x)$x)) / sqrt(p/(n-1)))
    }
    backvar <- mad1^2
    
    avgx <- t(t(x) - colMeans(x))
    dv <- svd(avgx)$d
    eigval_dat <- dv^2/(n-1)

    ##pad with 0s
    eigval_dat <- c(eigval_dat, rep(0, p-length(eigval_dat)))
    eigval_sim <- eigval_dat

    if (icovest == 1) { #use soft 
        taub <- 0
        tauu <- .soft_covest(eigval_dat, backvar)$tau
        etau <- (tauu-taub) / 100
        ids <- rep(0, 100)
        ## tune to determine whether hard/soft more appropriate
        for (i in 1:100) {
            taus <- taub + (i-1)*etau
            eigval_temp <- eigval_dat - taus
            eigval_temp[eigval_temp < backvar] <- backvar
            ids[i] <- eigval_temp[1] / sum(eigval_temp)
        }
        tau <- taub + (which.max(ids)-1)*etau
        eigval_sim <- eigval_dat - tau
        eigval_sim[eigval_sim < backvar] <- backvar

    } else if (icovest == 2) { #use sample eigenvalues
        eigval_sim[eigval_sim < 0] <- 0
        
    } else if (icovest == 3) { #use hard thresholding
        eigval_sim[eigval_dat < backvar] <- backvar
    } else {
        stop("covest must be 1, 2 or 3")
    }

    list(eigval_dat = eigval_dat,
         backvar = backvar,
         eigval_sim = eigval_sim)
}



## helper function for computing soft thresholding estimator
.soft_covest <- function(vsampeigv, sig2b) {

    p <- length(vsampeigv)
    vtaucand <- vsampeigv - sig2b
    
    ##if all eigenvals > sig2b, just use sample eigenvals
    if (vtaucand[p] > 0) {
        return(list(veigvest = vsampeigv,
                    tau = 0))
    }  
    
    ##if not enough power, just use flat est as in Matlab impl
    if (sum(vsampeigv) <= p*sig2b) {
        return(list(veigvest = rep(sig2b, p),
                    tau = 0))
    }

    ##find threshold to preserve power
    which <- which(vtaucand <= 0)
    icut <- which[1] - 1
    powertail <- sum(vsampeigv[(icut+1):p])
    power2shift <- sig2b*(p-icut) - powertail
    
    vi <- c(1:icut)
    vcumtaucand <- sort(cumsum(sort(vtaucand[vi])), decreasing=TRUE)
    
    vpowershifted <- (vi-1)*vtaucand[vi] + vcumtaucand

    flag <- (vpowershifted < power2shift)
    if (sum(flag) == 0) {
        ## means decreasing everything still not enough
        itau <- 0
    } else {
        ## means for some index, decrease is sufficient
        itau <- which(flag)[1]
    }

    if (itau == 1) {
        powerprop <- power2shift/vpowershifted[1] #originally no [1] idx, PKK
        tau <- powerprop*vtaucand[1]
    } else if (itau == 0) {
        powerprop <- power2shift/vpowershifted[icut] 
        tau <- powerprop*vtaucand[icut] 
    } else {
        powerprop <- (power2shift-vpowershifted[itau]) /
            (vpowershifted[itau-1]-vpowershifted[itau])
        tau <- vtaucand[itau] + powerprop*(vtaucand[itau-1] - vtaucand[itau])
    }

    
    veigvest <- vsampeigv - tau 
    flag <- (veigvest > sig2b) 
    veigvest <- flag*veigvest + (1-flag)*(sig2b*rep(1, p))

    ##return eigenvalue estimate and soft threshold parameter, tau
    list(veigvest = veigvest,
         tau = tau)
}










