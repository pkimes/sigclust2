#' Statistical Significance of Clustering (sigclust)
#' 
#' Re-implementation of the Monte Carlo simulation based
#' significance testing procedure described in Liu et al. (2008).
#' 
#' @param x a dataset with n rows and p columns, with observations in rows and
#'        features in columns.
#' @param n_sim a numeric value specifying the number of simulations for Monte Carlo 
#'        testing. (default = 100)
#' @param n_start a numeric value specifying the number of random starts to be used
#'        for k-means clustering passed to \code{kmeans}. (default = 1)
#' @param icovest an integer (1, 2 or 3) specifying the null covariance
#'        estimation method to be used. See \code{\link{null_eigval}} for more
#'        details. (default = 1)
#' @param bkgd_pca a logical value specifying whether to use scaled PCA scores
#'        over raw data to estimate the background noise. When FALSE, raw estimate
#'        is used; when TRUE, minimum of PCA and raw estimates is used.
#'        (default = FALSE)
#' @param labels a n-vector of 1s and 2s specifying cluster labels for testing
#'        instead of using \code{kmeans}. (default = NULL)
#' 
#' @return
#' The function returns a \code{sigclust} S3-object containing the 
#' resulting p-values. The \code{sigclust} object has following attributes:
#' \itemize{
#' \item{\code{in_mat}}: {the original data matrix passed to the constructor}
#' }
#' 
#' @references
#' \itemize{
#'     \item Liu Y., Hayes, D. N., Nobel, A., and Marron, J. S. (2008)
#'           Statistical significance of clustering for high-dimension, low-sample size data.
#'           Journal of the American Statistical Association.
#' }
#'
#' @export
#' @name sigclust
#' @aliases sigclust-constructor
#' @author Patrick Kimes
sigclust <- function(x, n_sim = 100, n_start = 1, icovest = 1,
                     bkgd_pca = FALSE, labels = NULL) {
    
    n <- nrow(x)
    p <- ncol(x)

    if (!is.matrix(x)) {
        stop("x must be a matrix; use as.matrix if necessary")
    }
    
    ## apply initial clustering
    dat <- .initcluster_sc2(x, n, p, n_start, labels)
    clust_dat <- dat$clust_dat
    ci_dat <- dat$ci_dat

    ##estimate null Gaussian
    x_null <- null_eigval(x, n, p, icovest, bkgd_pca)
    eigval_dat <- x_null$eigval_dat
    eigval_sim <- x_null$eigval_sim
    backvar <- x_null$backvar
    
    ## simulation indices
    ci_sim <- rep(0, n_sim)
    for (isim in 1:n_sim) {
        x_isim <- .simnull(eigval_sim, n, p)
        ci_sim[isim] <- .calc_simCI(x_isim, p, n_start)
    }

    ## report p-values
    p_norm <- pnorm(ci_dat, mean(ci_sim), sd(ci_sim))
    p_emp <- mean(ci_sim <= ci_dat)

    ## return sigclust structure
    structure(
        list(in_mat = x,
             in_args = list(n_sim = n_sim, n_start = n_start,
                 icovest = icovest, bkgd_pca = bkgd_pca),
             eigval_dat = eigval_dat,
             eigval_sim = eigval_sim,
             backvar = backvar,
             ci_dat = ci_dat,
             ci_sim = ci_sim,
             p_emp = p_emp,
             p_norm = p_norm),
        class = "sigclust")
}





## given null eigenvalues, simulate Gaussian dataset
.simnull <- function(eigval_sim, n, p) {
    simnorm <- matrix(rnorm(n*p, sd=sqrt(eigval_sim)), n, p, byrow=TRUE)
}

## calculate sum of squares
.sumsq <- function(x) { norm(sweep(x, 2, colMeans(x), "-"), "F")^2 }

## calculate 2-means cluster index (n x p matrices)
.calc2CI <- function(x1, x2) {
    if (is.matrix(x1) && is.matrix(x2) && ncol(x1) == ncol(x2)) {
        (.sumsq(x1) + .sumsq(x2)) / .sumsq(rbind(x1, x2))
    } else {
        stop(paste("x1, x2 must be matrices with same ncols",
                   "for 2CI calculation"))
    }      
}

## perform k-means clustering on the original data and 
## compute the corresponding cluster indices for each merge
.initcluster_sc2 <- function(x, n, p, n_start, labels) {
    if (is.null(labels)) {
        labels <- kmeans(x, 2, nstart=n_start)$cluster
    } else {
        ## check validity of labels passed to function
        if (length(labels) != n || any(sort(unique(labels)) != 1:2)) {
            stop("labels must be a n-vector of 1s and 2s")
        }
    }
    ci_dat <- .calc2CI(x[labels == 1, , drop=FALSE],
                       x[labels == 2, , drop=FALSE])
    list(ci_dat = ci_dat)
}

## given null eigenvalues, simulate Gaussian dataset
.simnull <- function(eigval_sim, n, p) {
    simnorm <- matrix(rnorm(n*p, sd=sqrt(eigval_sim)), n, p, byrow=TRUE)
}

## perform k-means clustering on a simulated dataset and
## compute the correspond cluster indices for only the final merge
.calc_simCI <- function(x, p, n_start) { 
    ## obtain clustering solution
    clust_isim <- kmeans(x, 2, nstart=n_start)
    ci_isim <- .calc2CI(x[clust_isim$cluster == 1, , drop=FALSE],
                        x[clust_isim$cluster == 2, , drop=FALSE])
    ci_isim
}


