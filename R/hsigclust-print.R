#' @title print hsigclust object
#' 
#' @description 
#' quick output to screen when hsigclust object is evaluated
#' 
#' @details
#' some details
#' 
#' @name hsigclust-print
#' @export 
#' @author Patrick Kimes

.print.hsigclust <- function(hsigclust, ...) {
  cat("This is better than printing everything. Don't y'all agree? \n")
  print( paste("number of p-values < 0.05: ",
               colSums(hsigclust@mpvalnorm<0.05)) )
}
setMethod("print", signature(x="hsigclust"),
          function(x, arg="all", ...) {
            .print.hsigclust(x, arg, ...)
          })
