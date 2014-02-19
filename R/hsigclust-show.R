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
  cat("This is better than showing everything. \n")
  print( paste("number of p-values < 0.05: ",
               colSums(hsigclust@mpvalnorm<0.05)) )
}
