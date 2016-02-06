#' provide summary of shc object
#'
#' Provides information about \code{shc} object similar to call to \code{show}
#'
#' @param object a \code{shc} object
#' @param ... other parameters to be used by the function
#' 
#' @export
#' @method summary shc
#' @name summary-shc
#' @author Patrick Kimes
summary.shc <- function(object, ...) {
    cat("shc object created using shc(..)\n\n")
    cat(paste0("clustered using:",
               "\n    linkage = ",
               object$in_args$linkage,
               "\n    dissimilarity = ",
               object$in_args$metric,
               "\n"))
    
    cat(paste0("SHC applied with:",
               "\n    n_sim = ",
               object$in_args$n_sim,
               "\n    icovest = ",
               object$in_args$icovest,
               "\n    ci = ",
               paste0(object$in_args$ci, collapse=", "),
               "\n    null_alg = ",
               paste0(object$in_args$null_alg, collapse=", "),
               "\n"))
    
    cat(paste0("FWER control:",
               "\n    ", 
               (object$in_args$alpha < 1),
               if (object$in_args$alpha < 1)
                   paste0(" [w/ alpha = ", 
                          object$in_args$alpha,
                          ", determined using CI# ",
                          object$in_args$ci_idx,
                          "]"),
               "\n    min obs for testing = ",
               object$in_args$n_min,
               "\n"))

    cat(paste0("# tests w/ p_norm < 0.05: ",
               colSums(object$p_norm < 0.05),
               ".\n\n") )
}
