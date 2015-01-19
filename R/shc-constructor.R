#' statistical Significance for Hierarchical Clustering (SHC) 
#' algorithm
#'
#' implements the Monte Carlo simulation based
#' significance testing procedure described in Kimes et al. (2014+).
#' Statistical significance is evaluated at each node along the tree
#' starting from the root using a Gaussian null hypothesis test. 
#' A corresponding family-wise error rate (FWER) controlling procedure is
#' provided.
#' 
#' @param x a dataset with n rows and p columns, with observations in rows
#' @param metric a string specifying the metric to be used in the hierarchical 
#'        clustering procedure. This must be a metric accepted by \code{dist}, 
#'        e.g. "euclidean," or "cor" (1 - Pearson correlation)
#'        (default = "euclidean")
#' @param linkage a string specifying the linkage to be used in the hierarchical 
#'        clustering procedure. This must be a linkage accepted by 
#'        \code{Rclusterpp.hclust} if \code{rcpp=TRUE}, e.g. "ward",
#'        or \code{stats::hclust}, if \code{rcpp=FALSE}, e.g. "ward.D2"
#'        (default = "ward.D2")
#' @param l an integer value specifying the power of the Minkowski distance, if 
#'        used (default = 2)
#' @param alpha a value between 0 and 1 specifying the desired level of the 
#'        test. If no FWER control is desired, simply set alpha to 1 (default = 1)
#' @param n_sim a numeric value specifying the number of simulations for Monte Carlo 
#'        testing (default = 100) 
#' @param n_min an integer specifying the minimum number of observations needed
#'        to calculate a p-value (default = 10)
#' @param icovest an integer between 1 and 3 specifying the null covariance
#'        estimation method to be used. See \code{\link{null_eigval}} for more
#'        details (default = 1)
#' @param bkgd_pca a logical value whether to use principal component scores when
#'        estimating background noise under the null (default = TRUE)
#' @param rcpp a logical value whether to use the \code{Rclusterpp} package
#'        (default = FALSE)
#' @param ci a string vector specifying the cluster indices to be used for 
#'        testing along the dendrogram. Currently, options include: "2CI", 
#'        "linkage". (default = "2CI")
#' @param null_alg a string vector specifying the clustering algorithms that 
#'        should be used to cluster the simulated null ditributions. Currently, options
#'        include: "2means" and "hclust". While "hclust" typically provides greater power
#'        for rotation invariant combinations of metric and linkage function, 
#'        if a non-rotation invariant metric, e.g. Pearson correlation, is used, it is
#'        recommended that the "2means" option is specified. Note, \code{null_alg} and
#'        \code{ci} must be of equal length. (default = "hclust")
#' @param ci_idx a numeric value between 1 and \code{length(ci)} 
#'        specifiying which CI to use for the FWER stopping rule.
#'        This only has an effect if \code{alpha} is specified to a non-default 
#'        value. (default = 1)
#' @param ci_emp a logical value specifying whether to use the empirical
#'        p-value from the CI based on \code{ci_idx} for the FWER stopping rule.
#'        As with \code{ci_idx}, this only has an effect if \code{alpha} is
#'        specified to a non-default value. (default = FALSE)
#' 
#' @return
#' The function returns a \code{shc} S3-object containing the 
#' resulting p-values. The \code{plot} method call will output a dendrogram
#' with the corresponding p-values placed at each merge. The \code{shc}
#' object has following attributes:
#' \itemize{
#' \item{\code{in_mat}}: {the original data matrix passed to the constructor}
#' \item{\code{in_args}}: {a list of the original parameters passed to the constructor}
#' \item{\code{eigval_dat}}: {a matrix with each row containing the sample eigenvalues for a
#'     subtree tested along the dendrogram}
#' \item{\code{eigval_sim}}: {a matrix with each row containing the estimated eigenvalues used to
#'     simulate null data at a subtree tested along the dendrogram}
#' \item{\code{backvar}}: {a vector containing the estimated background variances used for
#'     computing \code{eigval_sim}}
#' \item{\code{nd_type}}: {a vector of length n-1 taking values in "\code{n_small}",
#'     "\code{no_test}", "\code{tested}" specifying how each node along the
#'     dendrogram was handled by the iterative testing procedure}
#' \item{\code{ci_dat}}: {a matrix containing the cluster indices for the original
#'     data matrix passed to the constructor}
#' \item{\code{ci_sim}}: {a 3-dimensional array containing the simulated cluster
#'     indices at each subtree tested along the dendrogram}
#' \item{\code{p_emp}}: {a matrix containing the emprical p-values computed at each
#'     subtree tested along the dendrogram}
#' \item{\code{p_norm}}: {a matrix containing the Gaussian approximate p-values
#'     computed at each subtree tested along the dendrogram}
#' \item{\code{idx_hc}}: {a list of tuples containing the indices of clusters joined
#'     at each step of the hierarchical clustering procedure}
#' \item{\code{hc_dat}}: {a \code{hclust} object constructed from the original data matrix and
#'     arguments passed to the constructor}
#' }
#' 
#' @details
#' When possible, \code{Rclusterpp::Rclusterpp.hclust} should be used for
#' hierarchical clustering except in the case of clustering by Pearson
#' correlation (\code{dist="cor"}), for which we make use of
#' \code{WGCNA::cor} with the usual \code{stats::hclust}.
#'
#' For standard minimum variance Ward's linkage clustering, if \code{rcpp} is
#' \code{TRUE}, specify "ward" for \code{linkage}. However, if \code{rcpp} is
#' \code{FALSE}, then "ward.D2" should be specified. See \code{stats::hclust}
#' for details on changes to the function since R >= 3.0.4. 
#' 
#' The testing procedure will terminate when no nodes
#' meet the corresponding FWER control threshold specified by \code{alpha}
#' and \code{fwer}. Controlling the FWER  may also substantially speed up
#' the procedure by reducing the number of tests considered.
#'
#' If p-values are desired for all nodes, simply set \code{fwer} to \code{FALSE}
#' and \code{alpha} to 1. The FWER cutoffs may still be computed using
#' \code{fwer_cutoff} or by specifying \code{alpha} when calling \code{plot}.
#'
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 2), ncol = 2),
#'               matrix(rnorm(100, mean = -1), ncol = 2))
#' shc_cars <- shc(data, metric = "euclidean", linkage = "average")
#' tail(shc_cars$p_norm, 10)
#' 
#' @references
#' \itemize{
#'     \item Kimes, P. K., Hayes, D. N., Liu Y., and Marron, J. S. (2014)
#'           Statistical significance for hierarchical clustering.
#'           pre-print available.
#' }
#'
#' @export
#' @seealso \code{\link{plot-shc}} \code{\link{diagnostic}}
#' @import Rclusterpp WGCNA methods
#' @name shc
#' @aliases shc-constructor
#' @author Patrick Kimes
shc <- function(x, metric = "euclidean", linkage = "ward.D2", l = 2,
                alpha = 1, icovest = 1, bkgd_pca = TRUE, n_sim = 100,
                n_min = 10, rcpp = FALSE, ci = "2CI", null_alg = "hclust",
                ci_idx = 1, ci_emp = FALSE) {  
    
    n <- nrow(x)
    p <- ncol(x)

    n_ci <- length(ci)
    if (length(null_alg) != n_ci) {
        stop("ci and null_alg must be of same length")
    }

    for (ii in 1:n_ci) {
        if (ci[ii] == "linkage" && null_alg[ii] == "2means")
            stop("ci = 'linkage', null_alg = '2means' cannot be specified")
    }
        
    if (ci_idx > n_ci) {
        stop("invalid choice for ci_idx; ci_idx must be < length(ci)")
    }
    
    if (alpha > 1 || alpha < 0) {
        stop("invalid choice for alpha; alpha must be 0 < alpha < 1")
    }
    
    if (!is.matrix(x)) {
        stop("x must be a matrix; use as.matrix if necessary")
    }

    if (n_min < 3) {
        stop("n_min must be >= 3")
    }

    if (n_min > n) {
        stop("n_min must be <= n")
    }

    
    ##apply initial clustering
    x_clust <- .initcluster(x, n, p, metric, linkage, l, 
                           n_ci, ci, rcpp)
    ci_dat <- x_clust$ci_dat
    hc_dat <- x_clust$hc_dat
    idx_hc <- x_clust$idx_hc

    ##for plotting purposes, change heights of dendrogram
    if (linkage == "ward" & rcpp) {
        hc_dat$height <- sqrt(2*hc_dat$height)
    }
    
    ##p-values for all <=(n-1) tests
    p_emp <- matrix(2, nrow=n-1, ncol=n_ci)
    p_norm <- matrix(2, nrow=n-1, ncol=n_ci)
    colnames(p_emp) <- paste(null_alg, ci, sep="_")
    colnames(p_norm) <- paste(null_alg, ci, sep="_")
    
    ##null covariance parameters for all <=(n-1) tests
    eigval_dat <- matrix(-1, nrow=n-1, ncol=p)
    eigval_sim <- matrix(-1, nrow=n-1, ncol=p)
    backvar <- rep(-1, n-1)
    ci_sim <- array(-1, dim=c(n-1, n_sim, n_ci))
    
    ##determine parent nodes for all nodes
    pd_map <- .pd_map(hc_dat, n)
    
    ##compute Meinshausen cutoffs for significance at alpha
    cutoff <- fwer_cutoff(idx_hc, alpha)

    ##keep track of each node was tested
    nd_type <- rep("", n)
    nd_type[n] <- "sig"
    
    ##move through nodes of dendrogram in reverse order
    for (k in (n-1):1) { 

        ##indices for subtree
        idx_sub <- unlist(idx_hc[k, ])
        n_sub <- length(idx_sub)
        
        ##only calc p-values for branches w/ more than n_min
        if (n_sub < n_min) {
            nd_type[k] <- "n_small"
            next
        }

        ##if parent wasn't significant, skip
        ## - placed after n_min check on purpose
        if (nd_type[pd_map[k]] != "sig") {
            nd_type[k] <- "no_test"
            next 
        }

        ##estimate null Gaussian
        xk_null <- null_eigval(x[idx_sub, ], n_sub, p, icovest, bkgd_pca)

        ##prevent messages from looped application of clustering
        ## messages will be thrown once at initial clustering
        suppressMessages(
            ##simulate null datasets
            for (i in 1:n_sim) {
                xsim <- .simnull(xk_null$eigval_sim, n_sub, p)
                ci_sim[k, i, ] <- .calcCI_shc(xsim, p, metric, linkage, l, 
                                              n_ci, ci, null_alg, rcpp)
            }
        )
        
        ##compute p-values
        m_idx <- colMeans(as.matrix(ci_sim[k, , ]))
        s_idx <- apply(as.matrix(ci_sim[k, , ]), 2, sd)
        p_norm[k, ] <- pnorm(ci_dat[k, ], m_idx, s_idx)
        p_emp[k, ] <- colMeans(as.matrix(ci_sim[k, , ]) <= 
                                   matrix(ci_dat[k, ], nrow=n_sim, 
                                          ncol=n_ci, byrow=TRUE))

        ##flip p-values for linkage based testing
        p_norm[k, ci == "linkage"] <- 1-p_norm[k, ci == "linkage"]
        p_emp[k, ci == "linkage"] <- 1-p_emp[k, ci == "linkage"]

        ##keep everything
        eigval_dat[k, ] <- xk_null$eigval_dat
        eigval_sim[k, ] <- xk_null$eigval_sim
        backvar[k] <- xk_null$backvar
        
        ##update nd_type
        if (alpha < 1) {
            if (ci_emp) {
                nd_type[k] <- ifelse(p_emp[k, ci_idx] < cutoff[k],
                                     "sig", "not_sig")
            } else {
                nd_type[k] <- ifelse(p_norm[k, ci_idx] < cutoff[k],
                                     "sig", "not_sig")
            }
        } else {
            nd_type[k] <- "sig"
        }
        
    }
    
    structure(
        list(in_mat = x,
             in_args = list(metric = metric, linkage = linkage, alpha = alpha,
                 l = l, bkgd_pca = bkgd_pca, n_sim = n_sim,
                 n_min = n_min, icovest = icovest, ci = ci,
                 null_alg = null_alg, ci_idx = ci_idx, ci_emp = ci_emp),
             eigval_dat = eigval_dat,
             eigval_sim = eigval_sim,
             backvar = backvar,
             nd_type = nd_type[-n],
             ci_dat = ci_dat,
             ci_sim = ci_sim,
             p_emp = p_emp,
             p_norm = p_norm,
             idx_hc = idx_hc,
             hc_dat = hc_dat),
        class = "shc")
}



## #############################################################################
## #############################################################################
## helper functions

##identify parent node of each node in dendrogram
.pd_map <- function(hc, n) {
    ##determine parent branch node for all children nodes along dendrogram
    pd_pairs <- rbind(cbind(hc$merge[, 1], 1:(n-1)), 
                      cbind(hc$merge[, 2], 1:(n-1)))
    pd_map <- data.frame(pd_pairs[pd_pairs[, 1]>0, ])
    names(pd_map) <- c("dtr", "prt")
    pd_map <- pd_map[order(pd_map$dtr), 2] #the parent of each daughter
    pd_map <- c(pd_map, n) #add final node without a parent

    pd_map
}


##determine obs indices at each node of the denrogram
.idx_hc <- function(hc, n) {
    ##list array of cluster indices at each of the n-1 merges
    idx_hc <- array(list(), c(2*n-1, 2))
    idx_hc[1:n, 1] <- as.list(n:1)
    idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
    
    ##complete idx_hc
    for (k in 1:(n-1)) {
        idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
        idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
    }
    idx_hc <- idx_hc[-(1:n), ]

    idx_hc
}


##calculate sum of squares
.sumsq <- function(x) { norm(sweep(x, 2, colMeans(x), "-"), "F")^2 }


##calculate 2-means cluster index (n x p matrices)
.calc2CI <- function(x1, x2) {
    if (is.matrix(x1) && is.matrix(x2) && ncol(x1) == ncol(x2)) {
        (.sumsq(x1) + .sumsq(x2)) / .sumsq(rbind(x1, x2))
    } else {
        stop(paste("x1, x2 must be matrices with same ncols",
                   "for 2CI calculation"))
    }      
}
    

##parse clustering parameters to produce dendrogram
.cluster_shc <- function(x, metric, linkage, l, rcpp) {

    if (metric == "cor") {
        dmat <- 1 - WGCNA::cor(t(x))
        hc_dat <- hclust(as.dist(dmat), method=linkage)
    } else {
        if (rcpp) {
            hc_dat <- Rclusterpp.hclust(x, method=linkage, 
                                        distance=metric, p=l)
        } else {
            hc_dat <- hclust(dist(x, method=metric, p=l), method=linkage)
        }
    }
    hc_dat
}



##perform hierarchical clustering on the original data and 
## compute the corresponding cluster indices for each merge
.initcluster <- function(x, n, p, metric, linkage, l, 
                         n_ci, ci, rcpp) { 

    ##obtain clustering solution
    hc_dat <- .cluster_shc(x, metric, linkage, l, rcpp)

    ##list array of cluster indices at each of the n-1 nodes
    idx_hc <- .idx_hc(hc_dat, n)
    
    ##matrix containing cluster indices
    ci_dat <- matrix(-1, nrow=n-1, ncol=n_ci)

    ##calculate cluster index(ices) for merge k
    for (i_ci in 1:n_ci) {
        if (ci[i_ci] == "2CI") {
            for (k in 1:(n-1)) {
                ci_dat[k, i_ci] <- .calc2CI(x[idx_hc[[k, 1]], , drop=FALSE],
                                            x[idx_hc[[k, 2]], , drop=FALSE])
            }
        } else if (ci[i_ci] == "linkage") {
            ci_dat[, i_ci] <- hc_dat$height
        }
    }
    
    list(hc_dat = hc_dat, 
         idx_hc = idx_hc,
         ci_dat = ci_dat)
}


##given null eigenvalues, simulate Gaussian dataset
.simnull <- function(eigval_sim, n, p) {
    simnorm <- matrix(rnorm(n*p, sd=sqrt(eigval_sim)), n, p, byrow=TRUE)
}


##perform hierarchical clustering on a simulated dataset and
## compute the correspond cluster indices for only the final merge
.calcCI_shc <- function(x, p, metric, linkage, l, 
                        n_ci, ci, null_alg, rcpp) { 

    ##obtain clustering solution
    hc_isim <- .cluster_shc(x, metric, linkage, l, rcpp)
    split <- cutree(hc_isim, k=2)

    ##row vector containing cluster indices
    ci_isim <- matrix(-1, nrow=1, ncol=n_ci)

    for (i_ci in 1:n_ci) {
        if (ci[i_ci] == "2CI") {
            if (null_alg[i_ci] == "hclust") {
                ci_isim[i_ci] <- .calc2CI(x[split==1, , drop=FALSE],
                                          x[split==2, , drop=FALSE])        
            } else if (null_alg[i_ci] == "2means") {
                kmsol <- kmeans(x, centers=2)
                ci_isim[i_ci] <- kmsol$tot.withinss/kmsol$totss
            }
        } else if (ci[i_ci] == "linkage") {
            ci_isim[i_ci] <- hc_isim$height[nrow(x)-1]
        }
    }

    ci_isim
}









