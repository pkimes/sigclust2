#' null assumption diagnostic plots
#'
#' @name diagnostic-generic
#' @docType methods
#' @export
#' @keywords internal
diagnostic <- function(obj, ...) {
    UseMethod("diagnostic", obj)
}



#' plot diagnostics for shc object
#'
#' Provides visualizations to check the null Gaussian assumption at a specified
#' range of nodes along the dendrogram
#' 
#' @param obj a \code{shc} object
#' @param K an integer value specifying the number of nodes (starting from the
#'        root) for which to create diagnostic plots (default = 1)
#' @param fname a character string specifying the name of the output file,
#'        (default = "shc_diagnostic")
#' @param ... other parameters to be used by the function
#' 
#' @return
#' prints plots to \code{fname}.
#'
#' @export
#' @name diagnostic-shc
#' @author Patrick Kimes
diagnostic.shc <- function(obj, K = 1, fname = "shc_diagnostic", ...) {
    print("function not yet implemented. sorry.")
}
















{
    ##code taken from the original sigclust package
    ##
    ## .plot.sigclust <- function(sigclust, arg="all", ...) {
    ##   
    ##   raw.data <- sigclust@raw.data
    ##   veigval <- sigclust@veigval
    ##   simbackvar <- sigclust@simbackvar
    ##   vsimeigval <- sigclust@vsimeigval
    ##   icovest <- sigclust@icovest
    ##   nsim <- sigclust@nsim
  #   simcindex <- sigclust@simcindex
  #   pval <- sigclust@pval
  #   pvalnorm <- sigclust@pvalnorm
  #   xcindex <- sigclust@xcindex
  #   
  #   n <- dim(raw.data)[1]
  #   d <- dim(raw.data)[2]
  #   
  #   #Background Standard Deviation Diagnostic Plot
  #   
  #   overlay.x <- as.vector(raw.data)
  #   mean <- mean(overlay.x)
  #   sd <- sd(overlay.x)
  #   median <- median(overlay.x)
  #   mad <- mad(overlay.x)
  #   ntot <- length(overlay.x)
  #   maxnol <- 5000
  #   
  #   if(arg=="background"|arg=="all"){
  #     par(mfrow=c(1,1))
  #     par(mar=c(5,4,4,2)+0.1)
  #     nused <- maxnol
  #     denraw <- density(overlay.x)
  #     if(ntot>maxnol){
  #       overlay.x <- overlay.x[sample(c(1:ntot),maxnol)]
  #     }else{
  #       nused <- ntot
  #     }
  #     
  #     xmin <- min(denraw$x)
  #     xmax <- max(denraw$x)
  #     ymin <- min(denraw$y)
  #     ymax <- max(denraw$y)
  #     
  #     overlay.y <- ymin + (0.15+0.5*runif(nused))*(ymax - ymin)
  #     plot(denraw,xlim=range(overlay.x),ylim=range(denraw$y),
  #          col="blue",xlab="",main="",lwd=3,...)
  #     xgrid <- seq(xmin,xmax,by=0.0025*(xmax-xmin))
  #     normden <- dnorm(xgrid,mean=median,sd=sd)
  #     lines(xgrid,normden,col="red",lwd=3)
  #     points(overlay.x,overlay.y,col="green",pch=".")
  #     
  #     title("Distribution of All Pixel values combines")
  #     if(ntot>maxnol){
  #       text(xmin+0.47*(xmax-xmin),ymin+0.9*(ymax-ymin),
  #            paste("Overlay of",as.character(maxnol),"of",
  #                  as.character(ntot),"data points",sep=" "),cex=1.3)
  #     }else{
  #       text(xmin+0.47*(xmax-xmin),ymin+0.9*(ymax-ymin),
  #            paste("Overlay of", as.character(ntot),"data points",sep=" "),
  #            cex=1.3)
  #     }
  #     text(xmin+0.47*(xmax-xmin),ymin+0.8*(ymax-ymin),
  #          paste("Mean =",as.character(round(mean,3)),
  #                "  Median =",as.character(round(median,3)),sep=" "),cex=1.3)
  #     text(xmin+0.47*(xmax-xmin),ymin+0.7*(ymax-ymin),
  #          paste("s.d. =",as.character(round(sd,3)),
  #                "MAD =",as.character(round(mad,3)),sep=" "),cex=1.3)
  #     text(xmin+0.47*(xmax-xmin),ymin+0.6*(ymax-ymin),
  #          paste("Gaussian(",as.character(round(median,3)),",",
  #                as.character(round(mad,3)),") density",sep=""),
  #          col="red",cex=1.3)
  #     if(mad>sd){
  #       text(xmin+0.47*(xmax-xmin),ymin+0.55*(ymax-ymin),
  #            paste("Warning: MAD > s.d., SigClust can be anti-conservative",
  #                  sep=""),cex=1.3)
  #     }
  #   }
  #   if(arg=="qq"|arg=="all"){
  #     #QQ plot
  #     par(mfrow=c(1,1))  
  #     par(mar=c(5,4,4,2)+0.1)
  #     qqnorm <- qqnorm(as.vector(raw.data),plot.it=FALSE)
  #     if(ntot>maxnol){
  #       which <- sample(c(1:ntot),maxnol)
  #     }else{
  #       which <- c(1:ntot)
  #     }
  #     x <- sort(qqnorm$x[which])
  #     y <- sort(qqnorm$y[which])
  #     x25 <- x[which(x>qnorm(0.25))[1]]
  #     x50 <- x[which(x>qnorm(0.5))[1]]
  #     x75 <- x[which(x>qnorm(0.75))[1]]
  #     y25 <- y[which(x>qnorm(0.25))[1]]
  #     y50 <- y[which(x>qnorm(0.5))[1]]
  #     y75 <- y[which(x>qnorm(0.75))[1]]
  #     plot(x,y,col="red",xlab="Gaussian Q",ylab="Data Q",
  #          main="Robust Fit Gaussian Q-Q, All Pixel values",cex.lab=1.3,...)
  #     abline(0,1,col="green",lwd=2)
  #     xmin <- min(x)-0.05*(max(x)-min(x))
  #     xmax <- max(x)+0.05*(max(x)-min(x))
  #     ymin <- min(y)-0.05*(max(y)-min(y))
  #     ymax <- max(y)+0.05*(max(y)-min(y))
  #     text(xmin+0.3*(xmax-xmin),ymin+0.9*(ymax-ymin),
  #          paste("Mean =",as.character(round(mean,3)),sep=" "),cex=1.3)
  #     text(xmin+0.3*(xmax-xmin),ymin+0.8*(ymax-ymin),
  #          paste("sd =",as.character(round(sd,3)),sep=" "),cex=1.3)
  #     text(x25,y25,"+",cex=1.3)
  #     text(x25+0.7,y25,"0.25 quantile",cex=1.3)
  #     text(x50,y50,"+",cex=1.3)
  #     text(x50+0.7,y50,"0.5 quantile",cex=1.3)
  #     text(x75,y75,"+",cex=1.3)
  #     text(x75+0.7,y75,"0.75 quantile",cex=1.3)
  #   }
  #   if(arg=="diag"|arg=="all"){
  #     #Then make Covariance Estimation Diagnostic Plot
  #     
  #     
  #     ncut <- 100
  #     if(d>ncut){
  #       par(mfrow=c(2,2))
  #     }else{
  #       par(mfrow=c(1,2))
  #     }
  #     par(mar=c(2,3.7,2,1.7))
  #     veigvalpos <- veigval[which(veigval > 10^(-12))]
  #     dpos <- length(veigvalpos)
  #     xmin <- 0
  #     xmax <- d+1
  #     ymin <- min(veigval) - 0.05*(max(veigval)-min(veigval))
  #     ymax <- max(veigval) + 0.05*(max(veigval)-min(veigval))
  #     plot(c(1:d),vsimeigval,type="l",lty=2,lwd=3,col="red",
  #          xlim=c(xmin,xmax),ylim=c(ymin,ymax),
  #          xlab="Component #",ylab="Eigenvalue",...)
  #     points(c(1:d),veigval,col="black")
  #     title("Eigenvalues")
  #     lines(c(ncut+0.5,ncut+0.5),c(ymin,ymax),col="green")
  #     
  #     if(icovest!=2){
  #       lines(c(0,d+1),c(simbackvar,simbackvar),col="magenta")
  #       text(xmin+0.45*(xmax-xmin),ymin+0.9*(ymax-ymin),
  #            paste("Background variance = ",
  #                  as.character(round(simbackvar,3)),sep=""),
  #            col="magenta")
  #     }
  #     text(xmin+0.45*(xmax-xmin),ymin+0.8*(ymax-ymin),
  #          "Eigenvalues for simulation",col="red")
  #     if(mad>sd){
  #       text(xmin+0.45*(xmax-xmin),ymin+0.65*(ymax-ymin),
  #            "Warning: MAD > s.d.",col="magenta")
  #     }
  #     
  #     ymin <- min(log10(veigvalpos))
  #     - 0.05*(max(log10(veigvalpos)) - min(log10(veigvalpos)))
  #     ymax <- max(log10(veigvalpos))
  #     + 0.05*(max(log10(veigvalpos)) - min(log10(veigvalpos)))
  #     plot(c(1:d),log10(vsimeigval),type="l",lty=2,lwd=3,col="red",
  #          xlim=c(xmin,xmax),ylim=c(ymin,ymax),
  #          xlab="Component #",ylab="log10(Eigenvalue)",...)
  #     points(c(1:dpos),log10(veigvalpos),col="black")
  #     title("log10 Eigenvalues")
  #     lines(c(ncut+0.5,ncut+0.5),c(ymin,ymax),col="green")
  #     
  #     if(icovest!=2){
  #       lines(c(0,d+1),log10(c(simbackvar,simbackvar)),col="magenta")
  #       text(xmin+0.45*(xmax-xmin),ymin+0.9*(ymax-ymin),
  #            paste("log10 Background variance = ",
  #                  as.character(round(log10(simbackvar),3)),sep=""),
  #            col="magenta")
  #     }
  #     text(xmin+0.45*(xmax-xmin),ymin+0.8*(ymax-ymin),
  #          "Eigenvalues for simulation",col="red")
  #     if(mad>sd){
  #       text(xmin+0.45*(xmax-xmin),ymin+0.65*(ymax-ymin),
  #            "SigClust may be Anti-iConservative",col="magenta")
  #     }
  #     if(length(veigval)>=ncut){
  #       xmin <- 0
  #       xmax <- ncut+1
  #       ymin <- min(veigval[1:ncut])
  #       - 0.05*(max(veigval[1:ncut]) - min(veigval[1:ncut])) 
  #       ymax <- max(veigval[1:ncut])
  #       + 0.05*(max(veigval[1:ncut]) - min(veigval[1:ncut])) 
  #       plot(c(1:ncut),vsimeigval[1:ncut],type="l",lty=2,lwd=3,col="red",
  #            xlim=c(xmin,xmax),ylim=c(ymin,ymax),
  #            xlab="Component #",ylab="Eigenvalue",...)
  #       points(c(1:ncut),veigval[1:ncut],col="black")
  #       title("Zoomed in version of above")
  #       
  #       if(icovest!=2){
  #         lines(c(0,d+1),c(simbackvar,simbackvar),col="magenta")
  #         text(xmin+0.45*(xmax-xmin),ymin+0.9*(ymax-ymin),
  #              paste("Background variance = ",
  #                    as.character(round(simbackvar,3)),sep=""),
  #              col="magenta")
  #       }
  #       text(xmin+0.45*(xmax-xmin),ymin+0.8*(ymax-ymin),
  #            "Eigenvalues for simulation",col="red")
  #       
  #       nmax <- min(dpos,ncut)
  #       ymin <- min(log10(veigvalpos[1:nmax])) 
  #       - 0.05*(max(log10(veigvalpos[1:nmax])) - min(log10(veigvalpos[1:nmax])))
  #       ymax <- max(log10(veigvalpos[1:nmax])) 
  #       + 0.05*(max(log10(veigvalpos[1:nmax])) - min(log10(veigvalpos[1:nmax])))
  #       plot(c(1:nmax),log10(vsimeigval[1:nmax]),type="l",lty=2,lwd=3,
  #            col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
  #            xlab="Component #",ylab="log10(Eigenvalue)",...)
  #       points(c(1:nmax),log10(veigvalpos[1:nmax]),col="black")
  #       title("log10 Eigenvalues")
  #       
  #       if(icovest!=2){
  #         lines(c(0,ncut+1),log10(c(simbackvar,simbackvar)),col="magenta")
  #         text(xmin+0.30*(xmax-xmin),ymin+0.9*(ymax-ymin),
  #              paste("log10 Background variance = ",
  #                    as.character(round(log10(simbackvar),3)),sep=""),
  #              col="magenta")
  #       }
  #       text(xmin+0.30*(xmax-xmin),ymin+0.8*(ymax-ymin),
  #            "Eigenvalues for simulation",col="red")
  #     }
  #   }
  #   if(arg=="pvalue"|arg=="all"){
  #     
  #     # Make p Value plot
  #     
  #     par(mfrow=c(1,1))
  #     par(mar=c(5,4,4,2)+0.1)
  #     denpval <- density(simcindex)
  #     denrange <- quantile(denpval$y, probs=c(0,0.5, 0.75, 1))
  #     dy <- 0.1*(denrange[4]-denrange[1])
  #     
  #     mindex <- mean(simcindex)
  #     sindex <- sd(simcindex)
  #     
  #     xmin <- min(c(simcindex,xcindex))
  #     xmax <- max(c(simcindex,xcindex))
  #     dx <- 0.1*(xmax-xmin)
  #     xind <- seq(xmin,xmax,0.001)
  #     
  #     plot(denpval, xlim=c(xmin-dx,xmax+dx),col="red",
  #          xlab="Cluster Index", main="",lwd=2,...)
  #     title(main="SigClust Results")
  #     points(simcindex, runif(nsim, denrange[2], denrange[3]),
  #            col="blue",pch=".",cex=2)
  #     lines(c(xcindex, xcindex), c(denrange[1]-dy, denrange[4])+dy,
  #           col="green",lty=2,lwd=2)
  #     lines(xind,dnorm(xind,mean=mindex,sd=sindex),col="black",lty=3,lwd=2)
  #     legend(xmin+0.05*(xmax-xmin), denrange[4], paste("P-value=", pval),
  #            text.col="red",bty="n")
  #     legend(xmin+0.05*(xmax-xmin), denrange[3]+0.75*(denrange[4]-denrange[3]),
  #            paste("P-vNorm=", round(pvalnorm,3)), bty="n")
  #   }
  # }
  # 
  # 
  # setMethod("plot", signature(x="sigclust", y="missing"),
  #           function(x, y, arg="all", ...) {
  #             .plot.sigclust(x, arg, ...)
  #           })
}
