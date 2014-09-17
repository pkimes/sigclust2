#' provide summary of shc object
#'
#' Provides information about \code{shc} object similar to call to \code{show}
#'
#' @param object a \code{shc} object
#'
#' @export
#' @rdname summary-shc
#' @aliases summary summary,shc-method
#' @author Patrick Kimes
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
