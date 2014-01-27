#' @title show hsigclust object
#' 
#' @description 
#' quick output to screen when hsigclust object is evaluated.
#' I'm not exactly sure how this differs from the \code{print} method.
#' Currently under investigation, but defining both for safety.
#' 
#' @details
#' some details
#' 
#' @name hsigclust-show
#' @export 
#' @author Patrick Kimes

setMethod("show", signature(object="hsigclust"),
          function(object) {
            .show.hsigclust(object)
          })

.show.hsigclust <- function(hsigclust) {
  cat("This is better than SHOWing everything. Don't y'all agree? \n")
  print( paste("number of p-values < 0.05: ",
               colSums(hsigclust@mpvalnorm<0.05)) )
}
