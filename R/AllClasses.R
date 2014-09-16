#' significance for hierarchical clustering (shc) object
#'
#' S4 class encapsulating various elements of a SHC analysis.
#'
#' @slot in_mat the original data matrix passed to the constructor
#' @slot in_args a list of the original parameters passed to the constructor,
#'       see \code{shc} for list of parameters
#' @slot eigval_dat a matrix containing the sample eigenvalues for each subtree
#'       tested along the dendrogram
#' @slot eigval_sim a matrix containing the estimated eigenvalues used to
#'       simulate null data at each subtree tested along the dendrogram
#' @slot backvar a vector containing the estimated background variances used for
#'       computing \code{eigval_sim}
#' @slot nd_type a vector of length n-1 taking values in "\code{n_small}",
#'       "\code{no_test}", "\code{tested}" specifying how each node along the
#'       dendrogram was handled by the iterative testing procedure
#' @slot ci_dat a matrix containing the cluster indices for the original
#'       data matrix passed to the constructor
#' @slot ci_sim a 3-dimensional array containing the simulated cluster
#'       indices at each subtree tested along the dendrogram
#' @slot pval_emp a matrix containing the emprical p-values computed at each
#'       subtree tested along the dendrogram
#' @slot pval_norm a matrix containing the Gaussian approximate p-values
#'       computed at each subtree tested along the dendrogram
#' @slot idx_hc a list of tuples containing the indices of clusters joined
#'       at each step of the hierarchical clustering procedure
#' @slot hc_dat a \code{hclust} object constructed from the original data matrix and
#'       arguments passed to the constructor
#' 
#' @param obj \code{shc} object
#'
#' @exportClass shc
#' @name shc-class
#' @aliases shc-class
#' @rdname shc-class
#' @author Patrick Kimes

setOldClass("shc")
setClass("shc",
         slots=list(in_mat = "matrix",
                    in_args = "list",
                    eigval_dat = "matrix",
                    eigval_sim = "matrix",
                    backvar = "vector",
                    nd_type = "vector",
                    ci_sim = "array",
                    ci_dat = "matrix",
                    p_emp = "matrix",
                    p_norm = "matrix",
                    idx_hc = "array",
                    hc_dat = "hclust")
         )
