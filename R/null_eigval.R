#' Eigenvalue estimation for shc/sigclust testing procedures
#'
#' function to compute the eigenvalues of the null Gaussian distribution for
#' 
#' @param x a matrix of size n by p containing the original data
#' @param n an integer number of samples
#' @param p an integer number of features/covariates
#' @param icovest an integer between 1 and 3 corresponding to the covariance
#'        estimation procedure to use. See details for more
#'        information on the possible estimation procedures (default = 1)
#' @param bkgd_pca a logical value specifying whether to use scaled PCA scores
#'        or raw data to estimate the background noise (default = TRUE)
#' 
#' @return
#' The function returns a \code{hsigclust} object containing the 
#' resulting p-values. The print method call will output a dendrogram
#' with the corresponding p-values placed at each merge. 
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
#' Huang, H., Liu, Y., Yuan, M., and Marron, J. S. (2014).
#'  Statistical Significance of Clustering using Soft Thresholding.
#'  Journal of Computational and Graphical Statistics, preprint. \n
#' Liu, Y., Hayes, D. N., Nobel, A. B., and Marron, J. S. (2008).
#'  Statistical Significance of Clustering for High-Dimension, Low-Sample Size Data.
#'  Journal of the American Statistical Association, 103(483):1281–1293.
#' 
#' @export null_eigval
#' @name null_eigval
#' @author Patrick Kimes

null_eigval <- function(x, n, p, icovest = 1, bkgd_pca = TRUE) {

    if (nrow(x) != n | ncol(x) != p)
        stop("Wrong size of matrix x!")
    
    ## compute background based on PCA scores or raw data
    if (bkgd_pca) {
        mad1 <- mad(as.matrix(prcomp(x)$x)) / sqrt(p/(n-1))
    } else {
        mad1 <- mad(as.matrix(x))
    }
    simbackvar <- mad1^2
    
    avgx <- t(t(x) - colMeans(x))
    dv <- svd(avgx)$d
    veigval <- dv^2/(n-1)
    vsimeigval <- veigval
    
    if (icovest == 1) { #use soft 
        taub <- 0
        tauu <- .soft_covest(veigval, simbackvar)$tau
        etau <- (tauu-taub)/100
        ids <- rep(0, 100)
        for(i in 1:100){
            taus = taub + (i-1)*etau
            eigval.temp <- veigval - taus
            eigval.temp[eigval.temp<simbackvar] <- simbackvar
            ids[i] <- eigval.temp[1]/sum(eigval.temp)
        }
        tau <- taub + (which.max(ids)-1)*etau
        vsimeigval <- veigval - tau
        vsimeigval[vsimeigval<simbackvar] <- simbackvar

    } else if (icovest == 2) { #use sample eigenvalues
        vsimeigval[vsimeigval < 0] <- 0
        
    } else if (icovest == 3) { #use hard thresholding
        vsimeigval[veigval<simbackvar] <- simbackvar
    }

    list(veigval=veigval,
         simbackvar=simbackvar,
         vsimeigval=vsimeigval)
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
        itau <- 0
    } else {
        which <- which(flag > 0)
        itau <- which[1]
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










