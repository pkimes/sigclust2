.show.shc <- function(shc) {
    cat("shc object created using shc(..)\n\n")
    cat(paste0("clustered using:",
               "\n    linkage = ",
               shc@in_args$linkage,
               "\n    dissimilarity = ",
               shc@in_args$metric,
               "\n"))
    
    cat(paste0("SHC applied with:",
               "\n    n_sim = ",
               shc@in_args$n_sim,
               "\n    icovest = ",
               shc@in_args$icovest,
               "\n    ci = ",
               paste0(shc@in_args$ci, collapse=", "),
               "\n    ci_null = ",
               paste0(shc@in_args$ci_null, collapse=", "),
               "\n"))
    
    cat(paste0("FWER control:",
               "\n    ", 
               (shc@in_args$alpha < 1),
               if (shc@in_args$alpha < 1)
                   paste0(" [w/ alpha = ", 
                          shc@in_args$alpha,
                          ", determined using CI# ",
                          shc@in_args$ci_idx,
                          "]"),
               "\n    min obs for testing = ",
               shc@in_args$n_min,
               "\n"))
}



#' @describeIn shc show brief summary of \code{shc} object
#' @aliases show,shc-method
setMethod("show", signature(object="shc"),
          function(object) {
              .show.shc(object)
          })

