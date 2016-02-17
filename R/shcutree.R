#' Cut SHC Tree by Significance
#'
#' Returns the cluster labels based on cutting the tree at the specified
#' \code{alpha} level 
#' 
#' @param obj \code{shc} object
#' @param alpha a double between 0 and 1 specifying the FWER controlled significance
#'        cutoff (default = 0.05)
#' @param ci_idx a numeric value between 1 and \code{length(ci)} 
#'        specifiying which CI to use for the FWER stopping rule.
#'        This only has an effect if \code{alpha} < 1. (default = 1)
#' @param ci_emp a logical value specifying whether to use the empirical
#'        p-value from the CI based on \code{ci_idx} for the FWER stopping rule.
#'        As with \code{ci_idx} this only has an effect if \code{alpha} < 1.
#'        (default = FALSE)
#'
#' @return
#' a vector of length n specifying cluster labels based on 
#' 
#' @name shcutree
#' @export
#' @author Patrick Kimes
shcutree <- function(obj, alpha = 0.05, ci_idx = 1, ci_emp = FALSE) {
    
    ## check validity of alpha
    if (alpha > 1 || alpha < 0) {
        stop("invalid choice for alpha; alpha must be 0 < alpha < 1")
    }

    ## if method originally implemented with fwer control
    ## must use FWER and same ci_emp, ci_idx
    if (obj$in_args$alpha < 1) {
        if (alpha > obj$in_args$alpha) {
            warning(paste0("shc constructed using smaller alpha ",
                           "than specified; using alpha = ", obj$in_args$alpha))
        }
        if (obj$in_args$ci_idx != ci_idx) {
            warning("shc constructed using alpha < 1, using ci_idx from obj$in_args")
        }
        if (obj$in_args$ci_emp != ci_emp) {
            warning("shc constructed using alpha < 1, using ci_emp from obj$in_args")
        }
        alpha <- min(alpha, obj$in_args$alpha)
        ci_idx <- obj$in_args$ci_idx
        ci_emp <- obj$in_args$ci_emp
    }

    ## for easier calling, subset and reverse order
    if (ci_emp) {
        p_use <- rev(obj$p_emp[, ci_idx])
    } else {
        p_use <- rev(obj$p_norm[, ci_idx])
    }
    
    n <- length(p_use)  + 1

    ## determine significant clusters
    cutoff <- rev(fwer_cutoff(obj, alpha))
    pd_map <- .pd_map(obj$hc_dat, n)
    idx_hc <- obj$idx_hc[(n-1):1, ]
    
    nd_type <- rep("", n)
    nd_type[n] <- "sig"
    for (k in seq(n-1, by=-1)) {
        ## check if subtree is large enough
        if (length(unlist(idx_hc[k, ])) < obj$in_args$n_min) {
            nd_type[k] <- "n_small"
            next
        }
        ## check if parent was significant
        if (nd_type[pd_map[k]] != "sig") {
            nd_type[k] <- "no_test"
            next
        }
        ## compare against p-value cutoff
        if (alpha < 1) {
            nd_type[k] <- ifelse(p_use[k] < cutoff[k], "sig", "not_sig")
        }
    }

    ## determine which nodes have significant parent nodes
    ## but also aren't significant
    cl_idx <- which((pd_map %in% which(nd_type == "sig")) & (nd_type[-n] != "sig"))
    clusters <- apply(idx_hc[cl_idx, , drop=FALSE], 1, unlist)

    ## verify that the clusters contain all observations
    if (length(unique(unlist(clusters))) != n) {
        stop("length of cluster labels != n")
    }

    ## construct cluster label vector
    labs <- rep(-1, n)
    for (i in 1:length(clusters)) {
        labs[clusters[[i]]] <- i
    }

    ## return labels
    labs
}
