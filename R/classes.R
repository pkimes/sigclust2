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
#causes error with roxygen2 if setOldClass() function is at top



#S3 objects are a lot easier to work with... maybe in a later project.
# for the sake of knowledge, below is how one (I) would implement hsigclust
# as an S3 class.
# References:
#  http://www.pitt.edu/~njc23/Lecture4.pdf
hsigclust <- function(raw.data, meigval, msimeigval, vsimbackvar, icovest,
                      nsim, asimcindex, mpval, mpvalnorm, xmcindex,
                      clusterList, hc) {
  out <- list(raw.data=raw.data, 
              meigval=meigval, 
              msimeigval=msimeigval, 
              vsimbackvar=vsimbackvar, 
              icovest=icovest,
              nsim=nsim,
              asimcindex=asimcindex,
              mpval=mpval,
              mpvalnorm=mpvalnorm,
              xmcindex=xmcindex,
              clusterList=clusterList,
              hc=hc)
  class(out) <- "hsigclust"
  invisible(out)
}

