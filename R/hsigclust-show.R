#' @title show hsigclust object
#' 
#' @description Quick output to screen when hsigclust object is evaluated.
#' 
#' @details some details...
#' 
#' @name hsigclust-show
#' @export 
#' @author Patrick Kimes

setMethod("show", signature(object="hsigclust"),
          function(object) {
            .show.hsigclust(object)
          })

.show.hsigclust <- function(hsigclust) {
  cat("hsigclust object created using HSCtest()\n\n")
  cat(paste0("clustered using:",
             "\n    linkage = ",
             hsigclust@inparams$linkage,
             "\n    dissimilarity = ",
             hsigclust@inparams$metric,
             "\n"))
  
  cat(paste0("HSigClust applied with:",
             "\n    nsim = ",
             hsigclust@inparams$nsim,
             "\n    icovest = ",
             hsigclust@inparams$icovest,
             "\n    testCIs = ",
             paste0(hsigclust@inparams$testCIs, collapse=", "),
             "\n    testNulls = ",
             paste0(hsigclust@inparams$testNulls, collapse=", "),
             "\n"))
  
  cat(paste0("FWER control:",
             "\n    ", 
             hsigclust@inparams$alphaStop<1,
             if(hsigclust@inparams$alphaStop<1)
               paste0(" [w/ alpha = ", 
                      hsigclust@inparams$alphaStop,
                      ", determined using CI# ",
                      hsigclust@inparams$cutoffCI,
                      "]"),
             "\n    min obs for testing = ",
             hsigclust@inparams$minObs,
             "\n"))
}
