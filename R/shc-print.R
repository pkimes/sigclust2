#' print description of shc object
#'
#' Provides information about \code{shc}
#'
#' @param x a \code{shc} object
#' @param ... other parameters to be used by the function
#' 
#' @export
#' @method print shc
#' @name print-shc
#' @author Patrick Kimes
print.shc <- function(x, ...) {
    cat("\n")
    cat("shc object created using shc(..)")
    cat("\n")
    cat("--------------------------------")
    cat("\n")
    cat(paste0("Clustering Parameters:",
               "\n    dissimilarity = ",
               ifelse(is.null(x$in_args$metric),
                      "custom", x$in_args$metric),
               "\n    linkage = ",
               x$in_args$linkage,
               "\n"))
    
    cat(paste0("Testing Parameters:",
               "\n    n_sim = ",
               x$in_args$n_sim,
               "\n    icovest = ",
               x$in_args$icovest,
               "\n    ci = ",
               paste0(x$in_args$ci, collapse=", "),
               "\n    null_alg = ",
               paste0(x$in_args$null_alg, collapse=", "),
               "\n    n_min = ",
               x$in_args$n_min,
               "\n"))
    
    cat(paste0("    FWER control = ",
               (x$in_args$alpha < 1),
               if (x$in_args$alpha < 1)
                   paste0("\n        alpha = ", x$in_args$alpha,
                          "\n        ci_idx = ", x$in_args$ci_idx),
               "\n\n"))
}
