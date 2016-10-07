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
#' a vector of length n specifying cluster labels based on specified FWER
#' alpha threshold.
#' 
#' @name shcutree
#' @export
#' @author Patrick Kimes
shcutree <- function(obj, alpha = 0.05, ci_idx = 1, ci_emp = FALSE) {
    
    ## check validity of alpha
    if (alpha > 1 || alpha < 0) {
        stop("invalid choice for alpha; alpha must be 0 < alpha < 1")
    }

    ## check validity of CI index
    if (!(ci_idx %in% 1:ncol(obj$p_emp))) {
        stop("invalid choice for ci_idx; ci_idx must be a valid column index from p_emp or p_norm")
    }

    ## if shc originally used with alpha must use FWER and same ci_emp, ci_idx
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
        p_use <- obj$p_emp[, ci_idx]
    } else {
        p_use <- obj$p_norm[, ci_idx]
    }
    
    n <- length(p_use)  + 1

    ## determine significant clusters
    cutoff <- fwer_cutoff(obj, alpha)
    pd_map <- .pd_map(obj$hc_dat, n)
    
    nd_type <- rep("", n-1)
    for (k in 1:(n-1)) {
        ## check if subtree is large enough
        if (length(unlist(obj$idx_hc[k, ])) < obj$in_args$n_min) {
            nd_type[k] <- "n_small"
            next
        }
        ## check if parent was significant
        if ((k > 1) && (nd_type[pd_map[k]] != "sig")) {
            nd_type[k] <- "no_test"
            next
        }
        ## compare against p-value cutoff
        if (alpha < 1) {
            nd_type[k] <- ifelse(p_use[k] < cutoff[k], "sig", "not_sig")
        }
    }

    ## determine non-sig nodes with sig parent nodes
    cl_idx <- which((pd_map %in% which(nd_type == "sig")) & (nd_type != "sig"))
    if (length(cl_idx) == 0) { cl_idx <- 1 }
    clusters <- lapply(apply(obj$idx_hc[cl_idx, , drop=FALSE], 1, c), unlist)

    ## verify that the clusters contain all observations
    if (length(unique(unlist(clusters))) != n) {
        stop("length of cluster labels != n")
    }

    ## construct cluster label vector
    labs <- rep(-1, n)
    for (i in 1:length(clusters)) {
        labs[clusters[[i]]] <- i
    }

    labs
}
