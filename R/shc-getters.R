
#' @describeIn shc get \code{in_mat} slot
#' @aliases in_mat,shc-method
setMethod("in_mat", "shc", function(obj) return(obj@in_mat))

#' @describeIn shc get \code{in_args} slot
#' @aliases in_args,shc-method
setMethod("in_args", "shc", function(obj) return(obj@in_args))

#' @describeIn shc get \code{eigval_dat} slot
#' @aliases eigval_dat,shc-method
setMethod("eigval_dat", "shc", function(obj) return(obj@eigval_dat))

#' @describeIn shc get \code{eigval_sim} slot
#' @aliases eigval_sim,shc-method
setMethod("eigval_sim", "shc", function(obj) return(obj@eigval_sim))

#' @describeIn shc get \code{backvar} slot
#' @aliases backvar,shc-method
setMethod("backvar", "shc", function(obj) return(obj@backvar))

#' @describeIn shc get \code{nd_type} slot
#' @aliases nd_type,shc-method
setMethod("nd_type", "shc", function(obj) return(obj@nd_type))

#' @describeIn shc get \code{ci_sim} slot
#' @aliases ci_sim,shc-method
setMethod("ci_sim", "shc", function(obj) return(obj@ci_sim))

#' @describeIn shc get \code{ci_dat} slot
#' @aliases ci_dat,shc-method
setMethod("ci_dat", "shc", function(obj) return(obj@ci_dat))

#' @describeIn shc get \code{p_emp} slot
#' @aliases p_emp,shc-method
setMethod("p_emp", "shc", function(obj) return(obj@p_emp))

#' @describeIn shc get \code{p_norm} slot
#' @aliases p_norm,shc-method
setMethod("p_norm", "shc", function(obj) return(obj@p_norm))

#' @describeIn shc get \code{idx_hc} slot
#' @aliases idx_hc,shc-method
setMethod("idx_hc", "shc", function(obj) return(obj@idx_hc))

#' @describeIn shc get \code{hc_dat} slot
#' @aliases hc_dat,shc-method
setMethod("hc_dat", "shc", function(obj) return(obj@hc_dat))


#' @describeIn shc return FWER cutoffs based on \code{idx_hc} slot
#' @aliases fwer_cutoff,shc,numeric-method
setMethod("fwer_cutoff", signature(obj="shc", alpha="numeric"),
          function(obj, alpha) { fwer_cutoff(obj@idx_hc, alpha) })

#' compute cutoffs based on hierachical clustering output from idx_hc
#' @keywords internal
#' @aliases fwer_cutoff,list,numeric-method
setMethod("fwer_cutoff", signature(obj="list", alpha="numeric"),
          function(obj, alpha) {
              alpha/(nrow(obj)+1) *
                  apply(obj, 1, function(x) { length(unlist(x)) })
          })
