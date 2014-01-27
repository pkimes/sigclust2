#' @title Validity check for the hsigclust class
#' 
#' @description
#' making sure the \code{hsigclust} object is correctly constructed
#' 
#' @details
#' some details
#' 
#' @references
#' http://www.bioconductor.org/help/course-materials/2010/AdvancedR/S4InBioconductor.pdf
#'  (slide 12)
#'  
#' @name setValidity
#' @author Patrick Kimes

#check validity of slot types
.validtype <- function(object) {
  if (!is.matrix(object@raw.data))
    return("'raw.data' is not a matrix")
  else if (!is.matrix(object@meigval))
    return("'meigval' is not a matrix")
  else if (!is.matrix(object@msimeigval))
    return("'msimeigval' is not a matrix")
  else if (!is.vector(object@vsimbackvar))
    return("'vsimbackvar' is not a vector")
  else if (!is.numeric(object@icovest))
    return("'icovest' is not a numeric")
  else if (!is.numeric(object@nsim)) #also check length==1
    return("'nsim' is not a numeric")
  else if (!is.array(object@asimcindex)) #check dims
    return("'asimcindex' is not an array")
  else if (!is.matrix(object@mpval))
    return("'mpval' is not a matrix")
  else if (!is.matrix(object@mpvalnorm))
    return("'mpvalnorm' is not a matrix")
  else if (!is.matrix(object@xmcindex))
    return("'xmcindex' is not a matrix")
  else if (!is.matrix(object@clusterList))
    return("'clusterList' is not a matrix")
  else if (!is.matrix(object@hc))
    return("'hc' is not a hclust object")
#   else
#     return(TRUE)
}

#check validity of slot dimensions (assuming data is correct size)
.validsize <- function(object) {
  if (!is.matrix(object@raw.data))
    return("'raw.data' is not a matrix")
  else if (!is.matrix(object@meigval))
    return("'meigval' is not a matrix")
  else if (!is.matrix(object@msimeigval))
    return("'msimeigval' is not a matrix")
  else if (!is.vector(object@vsimbackvar))
    return("'vsimbackvar' is not a vector")
  else if (!is.numeric(object@icovest))
    return("'icovest' is not a numeric")
  else if (!is.numeric(object@nsim)) #also check length==1
    return("'nsim' is not a numeric")
  else if (!is.array(object@asimcindex)) #check dims
    return("'asimcindex' is not an array")
  else if (!is.matrix(object@mpval))
    return("'mpval' is not a matrix")
  else if (!is.matrix(object@mpvalnorm))
    return("'mpvalnorm' is not a matrix")
  else if (!is.matrix(object@xmcindex))
    return("'xmcindex' is not a matrix")
  else if (!is.matrix(object@clusterList))
    return("'clusterList' is not a matrix")
  else if (!is.matrix(object@hc))
    return("'hc' is not a hclust object")
  #   else
  #     return(TRUE)  
}

#check validity of misc. properties
.validmisc <- function(object) {
  #1. check if hc has same data as slot raw.data
}

setValidity("hsigclust", 
            function(object) {
              .validtype(object)
              .validsize(object)
              .validmisc(object)
              return(TRUE) #since no errors caught above
            }
)

