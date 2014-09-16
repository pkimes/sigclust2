## ###################################################################
## ###################################################################
## constructor methods

#' @name shc
#' @export
#' @docType methods
#' @rdname shc-constructor
setGeneric("shc",
           valueClass="shc",
           function(x, ...) standardGeneric("shc"))



## ###################################################################
## ###################################################################
## misc methods

#' null assumption diagnostic plots
#'
#' Method for generating diagnostic plots for statistical significance
#' of clustering analysis
#' 
#' @name diagnostic
#' @export
#' @docType methods
#' @rdname diagnostic
setGeneric("diagnostic",
           function(obj, ...)  standardGeneric("diagnostic"))



## ###################################################################
## ###################################################################
## getter methods

#' in_mat
#' return \code{in_mat} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{in_mat}
#' @keywords internal
#' @rdname in_mat
setGeneric("in_mat", function(obj) standardGeneric("in_mat"))

#' in_args
#' return \code{in_args} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{in_args}
#' @keywords internal
#' @rdname in_args
setGeneric("in_args", function(obj) standardGeneric("in_args"))

#' eigval_dat
#' return \code{eigval_dat} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{eigval_dat}
#' @keywords internal
#' @rdname eigval_dat
setGeneric("eigval_dat", function(obj) standardGeneric("eigval_dat"))

#' eigval_sim
#' return \code{eigval_sim} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{eigval_sim}
#' @keywords internal
#' @rdname eigval_sim
setGeneric("eigval_sim", function(obj) standardGeneric("eigval_sim"))

#' backvar
#' return \code{backvar} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{backvar}
#' @keywords internal
#' @rdname backvar
setGeneric("backvar", function(obj) standardGeneric("backvar"))

#' nd_type
#' return \code{nd_type} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{nd_type}
#' @keywords internal
#' @rdname nd_type
setGeneric("nd_type", function(obj) standardGeneric("nd_type"))

#' ci_sim
#' return \code{ci_sim} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{ci_sim}
#' @keywords internal
#' @rdname ci_sim
setGeneric("ci_sim", function(obj) standardGeneric("ci_sim"))

#' ci_dat
#' return \code{ci_dat} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{ci_dat}
#' @keywords internal
#' @rdname ci_dat
setGeneric("ci_dat", function(obj) standardGeneric("ci_dat"))

#' p_emp
#' return \code{p_emp} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{p_emp}
#' @keywords internal
#' @rdname p_emp
setGeneric("p_emp", function(obj) standardGeneric("p_emp"))

#' p_norm
#' return \code{p_norm} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{p_norm}
#' @keywords internal
#' @rdname p_norm
setGeneric("p_norm", function(obj) standardGeneric("p_norm"))

#' idx_hc
#' return \code{idx_hc} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{idx_hc}
#' @keywords internal
#' @rdname idx_hc
setGeneric("idx_hc", function(obj) standardGeneric("idx_hc"))

#' hc_dat
#' return \code{hc_dat} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{hc_dat}
#' @keywords internal
#' @rdname hc_dat
setGeneric("hc_dat", function(obj) standardGeneric("hc_dat"))



## ###################################################################
## ###################################################################
## getter functions

#' fwer_cutoff
#' return Family-Wise Error Rate (FWER) cutoffs
#' @export
#' @docType methods
#' @param obj currently only for objects of class \code{shc-class}
#' @param alpha a numeric value between 0 and 1
#' @keywords internal
#' @rdname fwer_cutoff
setGeneric("fwer_cutoff", function(obj, alpha) standardGeneric("fwer_cutoff"))

