
#' @describeIn shc show brief summary of \code{shc} object
#' @aliases show,shc-method
setMethod("show", signature(object="shc"),
          function(object) {
            .show.shc(object)
          })

.show.shc <- function(shc) {
  cat("shc object created using shc(..)\n\n")
  cat(paste0("clustered using:",
             "\n    linkage = ",
             shc@in_args$linkage,
             "\n    dissimilarity = ",
             shc@in_args$metric,
             "\n"))
  
  cat(paste0("SHC applied with:",
             "\n    nsim = ",
             shc@in_args$nsim,
             "\n    icovest = ",
             shc@in_args$icovest,
             "\n    testCIs = ",
             paste0(shc@in_args$testCIs, collapse=", "),
             "\n    testNulls = ",
             paste0(shc@in_args$testNulls, collapse=", "),
             "\n"))
  
  cat(paste0("FWER control:",
             "\n    ", 
             (shc@in_args$alphaStop < 1),
             if (shc@in_args$alphaStop < 1)
               paste0(" [w/ alpha = ", 
                      shc@in_args$alphaStop,
                      ", determined using CI# ",
                      shc@in_args$cutoffCI,
                      "]"),
             "\n    min obs for testing = ",
             shc@in_args$minObs,
             "\n"))
}
