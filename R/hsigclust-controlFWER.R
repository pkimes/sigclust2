#' @title HSigClust FWER controlling procedure
#' 
#' @description
#' Calculate the sequence of p-value cutoffs for FWER control as described in
#' Meinshausen (2008)
#' 
#' @details 
#' some details
#' 
#' @references 
#' Meinshausen, N. (2008). Hierarchical testing of variable importance. Biometrika, 95(2), 265-278.
#'
#' @exportMethod controlFWER
#' @export
#' @name hsigclust-controlFWER 
#' @author Patrick Kimes


setGeneric("controlFWER", 
           function(hsc, alpha=0.05, ...) standardGeneric("controlFWER"))
setMethod("controlFWER", 
          signature(hsc="hsigclust"), 
          function(hsc, alpha, ...) .controlFWER.hsigclust(x, alpha, ...))

.controlFWER.hsigclust <- function(hsc, alpha=0.05, ...) {
  # function doesn't need entire vector. just cluster sizes,
  # e.g. attributes(hc[[1]])$members
  # 1. obtain cluster sizes at each step from hsigclust object hc(hsc)
  # 2. return the vector of 0.05/(clustersizes / #samples)
  # -- how do we incorporate this? - give it as user output?
  # -- have it be an option in original call to stop from top
  # -- other?
}


