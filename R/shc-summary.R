#' provide summary of shc object
#'
#' Provides information about \code{shc} object similar to call to \code{print}
#'
#' @param object a \code{shc} object
#' @param ... other parameters to be used by the function
#' 
#' @export
#' @method summary shc
#' @name summary-shc
#' @author Patrick Kimes
summary.shc <- function(object, ...) {
    cat("\n")
    cat("shc object created using shc(..)")
    cat("\n")
    cat("--------------------------------")
    cat("\n")
    cat(paste0("Clustering Parameters:",
               "\n    dissimilarity = ",
               ifelse(is.null(object$in_args$metric),
                      "custom", object$in_args$metric),
               "\n    linkage = ",
               object$in_args$linkage,
               "\n"))
    
    cat(paste0("Testing Parameters:",
               "\n    n_sim = ",
               object$in_args$n_sim,
               "\n    icovest = ",
               object$in_args$icovest,
               "\n    ci = ",
               paste0(object$in_args$ci, collapse=", "),
               "\n    null_alg = ",
               paste0(object$in_args$null_alg, collapse=", "),
               "\n    n_min = ",
               object$in_args$n_min,
               "\n"))
    
    cat(paste0("    FWER control = ",
               (object$in_args$alpha < 1),
               if (object$in_args$alpha < 1)
                   paste0("\n        alpha = ", object$in_args$alpha,
                          "\n        ci_idx = ", object$in_args$ci_idx),
               "\n\n"))
}
