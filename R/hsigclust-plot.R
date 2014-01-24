#' @title plot hsigclust object
#' 
#' @description 
#' visualize results of HSigClust analysis as an annotated dendrogram
#' object using the ggdendro and ggplot2 packages
#' 
#' @details
#' some details
#' 
#' @import ggplot2 ggdendro
#' @name hsigclust-plot
#' @export
#' @author Patrick Kimes

.plot.hsigclust <- function(hsigclust, arg="all", ...) {
  #first plot plain dendrogram
  plot(hsigclust@hc, xlab="", labels=FALSE, )
  ntest <- nrow(hsigclust@mpvalnorm)
  #below taken from pvclust code...
  usr <- par("usr")
  xwd <- usr[2] - usr[1]
  ywd <- usr[4] - usr[3]
  cin <- par()$cin
  for (k in ntest:1) {
    if (hsigclust@mpvalnorm[k, 1] < 0.05) {
      mi <- unlist(hsigclust@clusterList[k, ])
      ma <- match(mi, hsigclust@hc$order)
      
      xl <- min(ma)
      xr <- max(ma)
      mx <- xwd / 3 / ntest
      
      yt <- hsigclust@hc$height[k]
      yb <- usr[3]
      my <- ywd / 100
      
      rect(xl - mx, yt - my, xr + mx, yt + my, border=2, shade=NULL)
    } 
  } 
  {
  ########### from pvclust
  col <- c(2, 3, 8) 
  print.num <- TRUE
  float <- 0.01 
  cex <- .8
  font <- NULL
  nlabels <- 5
  plotidx <- tail(1:nrow(hsigclust@mpval), nlabels)
  
  axes <- hc2axes(hsigclust@hc)
  usr  <- par()$usr
  wid <- usr[4] - usr[3]
  pv <- format(signif(hsigclust@mpval[plotidx, 1], 2), scientific=TRUE)
  pvn <- format(signif(hsigclust@mpvalnorm[plotidx, 1], 2), scientific=TRUE)
  rn <- as.character(plotidx)
  pv[length(pv)] <- paste0("pval: \n", pv[length(pv)])
  pvn[length(pvn)] <- paste0("pvalnorm: \n", pvn[length(pvn)])
  rn[length(rn)] <- paste0("edge #: \n", rn[length(rn)])
  a <- text(x=axes[plotidx,1], y=axes[plotidx,2] + float * wid, pv,
            col=col[1], pos=2, offset=.3, cex=cex, font=font)
  a <- text(x=axes[plotidx,1], y=axes[plotidx,2] + float * wid, pvn,
            col=col[2], pos=4, offset=.3, cex=cex, font=font)
  if(print.num)
    a <- text(x=axes[plotidx,1], y=axes[plotidx,2], rn,
              col=col[3], pos=1, offset=.3, cex=cex, font=font)
  } 
}

setMethod("plot", signature(x="hsigclust", y="missing"),
          function(x, y, arg="all", ...) {
            .plot.hsigclust(x, arg, ...)
          })




# DIRECTLY FROM pvclust, will remove with update 
# to ggdendro plotting
hc2axes <- function(x)
{
  A <- x$merge # (n,n-1) matrix
  n <- nrow(A) + 1
  x.axis <- c()
  y.axis <- x$height
  
  x.tmp  <- rep(0,2)
  zz     <- match(1:length(x$order),x$order)
  
  for(i in 1:(n-1)) {
    ai <- A[i,1]
    
    if(ai < 0)
      x.tmp[1] <- zz[-ai]
    else
      x.tmp[1] <- x.axis[ai]
    
    ai <- A[i,2]
    
    if(ai < 0)
      x.tmp[2] <- zz[-ai]
    else
      x.tmp[2] <- x.axis[ai]
    
    x.axis[i] <- mean(x.tmp)
  }
  
  return(data.frame(x.axis=x.axis,y.axis=y.axis))
}

