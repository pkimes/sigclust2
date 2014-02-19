#' @title hsigclust class 
#' 
#' @name hsigclust
#' @aliases hsigclust-class
#' @rdname hsigclust
#' 
#' @description The hsigclust class is used to store output of the HSigClust 
#'              testing procedure implemented in \code{HSCtest}. 
#' 
#' @details some details...
#' 
#' @exportClass hsigclust
#' @author Patrick Kimes

setOldClass("hclust")
setClass("hsigclust",
         slots=list(raw.data = "matrix",
                    meigval = "matrix",
                    msimeigval = "matrix",
                    vsimbackvar = "vector",
                    icovest = "numeric",
                    nsim = "numeric",
                    asimcindex = "array",
                    mpval = "matrix",
                    mpvalnorm = "matrix",
                    xmcindex = "matrix",
                    clusterList = "matrix",
                    hc = "hclust")
         )
