#' @title hsigclust class 
#' 
#' @description 
#' the hsigclust class used to store output of a HSigClust analysis 
#' performed using the HSCtest() function
#' 
#' @details
#' some details
#' 
#' @name hsigclust
#' @rdname hsigclust
#' @exportClass hsigclust
#' @aliases hsigclust-class
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
