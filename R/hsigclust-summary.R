#' @title summary for hsigclust object
#' 
#' @description 
#' details on the results of the HSigClust testing procedure 
#' 
#' @details
#' some details
#' 
#' @name hsigclust-summary
#' @export 
#' @author Patrick Kimes

.summary.hsigclust <- function(hsigclust, arg="all", ...) {
  print("This is a horrible summary. Do something with it.")
}
setMethod("summary", signature(object="hsigclust"),
          function(object, ...) {
            .summary.hsigclust(object, arg="all", ...)
          })

