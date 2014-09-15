#' @describeIn shc show brief summary of shc object
#' @aliases summary,shc-method

setMethod("summary", signature(object="shc"),
          function(object, ...) {
            .summary.shc(object, ...)
          })

.summary.shc <- function(shc, ...) {
  show(shc)
  cat("\n")
  cat(paste0("number of p-values < 0.05: ",
            colSums(shc@p_norm < 0.05),
            ".\n\n") )
}
