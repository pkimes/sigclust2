#' Run Statistical Significance of Hierarchical Clustering algorithm
#'
#' @name HSCtest
#'
#' @description \code{HSCtest} first computes the desired dendrogram, and quantifies the 
#'              statistical significance of clustering at each branching along the tree using
#'              a hypothesis testing procedure. The code is written so that various cluster 
#'              inidices can be easily introduced by writting new CI function and adding 
#'              to .initcluster(), .simcluster() and incrementing "nCIs." 
#'              -- more elegant solution?
#' 
#' @param x a dataset with p rows and n columns, with observations in columns.
#' @param metric a string specifying the metric to be used in the hierarchical clustering procedure.
#'        This must be a metric accepted by \code{dist}, e.g. "euclidean." If squared Euclidean distance
#'        (or the square of any other metric) is desired, set the \code{square} parameter to \code{TRUE}.
#' @param linkage a string specifying the linkage to be used in the hierarchical clustering procedure.
#'        This must be a linkage accepted by \code{hclust}, e.g. "ward."
#' @param alpha a value between 0 and 1 specifying the desired level of the test. If no FWER control is
#'        desired, simply set alpha to 1. The default is 0.05.
#' @param square a logical specifying whether to square the dissimilarity matrix produced by specified 
#'        \code{metric}. This is necessary, for example, in order to implement Ward's minimum variance
#'        method with squared Euclidean metric. (Honestly can't think of any other situations when 
#'        you'd want to square the diss. matrix.) Default is FALSE. 
#' @param l an integer value specifying the power of the Minkowski distance, if used, default is 2.
#' @param nsim a numeric value specifying the number of simulations for SigClust testing.
#'        The default is to run 100 simulations at each merge. 
#' @param minObs an integer specifying the minimum number of observations needed to calculate a p-value,
#'        default is 10.
#' @param icovest a numeric value specifying the covariance estimation method: 
#'        1. Use a soft threshold method as constrained MLE (default); 
#'        2. Use sample covariance estimate (recommended when diagnostics fail); 
#'        3. Use original background noise thresholded estimate (from Liu et al., (2008)) 
#'        ("hard thresholding") as described in the \code{sigclust} package documentation.
#' 
#' @return The function returns a \code{hsigclust} object containing the resulting p-values.
#'        The print method call will output a dendrogram with the corresponding p-values placed 
#'        at each merge. 
#' 
#' 
#' @details The function extends the \code{sigclust} idea to the hierarchical setting by modifying the
#'          clustering procedure employed to compute the null distribution of cluster indices.
#' 
#' @import sigclust
#' @export HSCtest
#' @author Patrick Kimes

##need to run this since "hclust" isn't a formal S3 class
##ref: http://stackoverflow.com/questions/12636056/
#getClassDef("hclust")
#setOldClass("hclust")


#main function for performing HSigClust testing
HSCtest <- function(x, metric, linkage, alpha=0.05, square=FALSE, l=2, nsim=100, minObs=10, icovest=1) {  

  #number of cluster indices
  nCIs <- 1

  #convert boolean to 1,2
  square <- square+1
  #check the dimension of x to match n and p
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  x <- as.matrix(x)
  xclust <- .initcluster(x, n, p, metric, linkage, square, l, nCIs)
  xmcindex <- xclust$mcindex
  
  #p-values for all <=(n-1) tests
  mpval <- matrix(2, nrow=n-1, ncol=nCIs)
  mpvalnorm <- matrix(2, nrow=n-1, ncol=nCIs)

  #null covariance parameters for all <=(n-1) tests
  meigval <- matrix(-1, nrow=n-1, ncol=p)
  msimeigval <- matrix(-1, nrow=n-1, ncol=p)
  vsimbackvar <- rep(-1, n-1)
  asimcindex <- array(-1, dim=c(n-1, nsim, nCIs))
  
  for (k in 1:(n-1)) {

    subxIdx <- unlist(xclust$clusterList[k, ])
    subn <- length(subxIdx)
    #only calc p-values for branches w/ more than minObs
    if (subn >= minObs) {
      xvareigen <- .vareigen(x[subxIdx, ], subn, p, icovest)
      for (i in 1:nsim) {
        xsim <- .simnull(xvareigen$vsimeigval, subn, p, 
                         metric, linkage, square, l, nCIs)
        asimcindex[k, i, ] <- xsim$mcindex        
      }
      mindex <- colMeans(as.matrix(asimcindex[k, , ]))
      sindex <- apply(as.matrix(asimcindex[k, , ]), 2, sd)
      mpvalnorm[k, ] <- pnorm(xmcindex[k, ], mindex, sindex)
      mpval[k, ] <- colMeans(as.matrix(asimcindex[k, , ]) <= 
                               matrix(xmcindex[k, ], nrow=nsim, 
                                      ncol=nCIs, byrow=TRUE))
      meigval[k, ] <- xvareigen$veigval
      msimeigval[k, ] <- xvareigen$vsimeigval
      vsimbackvar[k] <- xvareigen$simbackvar
    }    
    
  }  
  
  return(new("hsigclust", 
             raw.data = x,
             meigval = meigval,
             msimeigval = msimeigval,
             vsimbackvar = vsimbackvar,
             icovest = icovest,
             nsim = nsim,
             asimcindex = asimcindex,
             mpval = mpval,
             mpvalnorm = mpvalnorm,
             xmcindex = xmcindex,
             clusterList = xclust$clusterList,
             hc = xclust$clusters))
  
}

#calculate 2-means cluster index
.sumsq <- function(x) { norm(sweep(x, 2, colMeans(x), "-"), "F")^2 }
.calc2CI <- function(x1, x2) {
  if (is.matrix(x1) && is.matrix(x2) && ncol(x1) == ncol(x2)) {
    ci <- (.sumsq(x1) + .sumsq(x2)) / .sumsq(rbind(x1, x2))
  } else {
    error(paste("x1, x2 must be matrices with same ncols",
                "for 2CI calculation"))
  }      
}

#perform hierarchical clustering on the original data and 
# compute the corresponding cluster indices for each merge
.initcluster <- function(x, n, p, metric, linkage, square, l, nCIs) { 

  #need to implement clustering algorithm
  dmatrix <- dist(x, method=metric, p=l)
  clusters <- hclust(dmatrix^square, method=linkage)

  #matrix containing cluster indices
  mcindex <- matrix(-1, nrow=n-1, ncol=nCIs)
  
  #list array of cluster indices at each of the n-1 merges
  clusterList <- array(list(), c(2*n-1, 2))
  clusterList[1:n, 1] <- as.list(n:1)
  clusterList[(n+1):(2*n-1), ] <- clusters$merge+n+(clusters$merge<0)
  for (k in 1:(n-1)) {
    clusterList[[n+k, 1]] <- unlist( clusterList[clusterList[[n+k, 1]], ] )
    clusterList[[n+k, 2]] <- unlist( clusterList[clusterList[[n+k, 2]], ] )
    #calculate cluster index(ices) for merge k
    mcindex[k, 1] <- .calc2CI(x[clusterList[[n+k, 1]], , drop=FALSE],
                           x[clusterList[[n+k, 2]], , drop=FALSE])
  }
  clusterList <- clusterList[-(1:n), ]
  
  return(list(clusters=clusters, 
              clusterList=clusterList,
              mcindex=mcindex))
}

#perform hierarchical clustering on a simulated dataset and
# compute the correspond cluster indices for only the final merge
.simcluster <- function(sim_x, p, metric, linkage, square, l, nCIs) { 
  #need to implement clustering algorithm
  dmatrix <- dist(sim_x, method=metric, p=l)
  clusters <- hclust(dmatrix^square, method=linkage)
  sim_split <- cutree(clusters, k=2)
  mcindex <- matrix(-1, nrow=1, ncol=nCIs)
  mcindex[1] <- .calc2CI(sim_x[sim_split==1, , drop=FALSE],
                         sim_x[sim_split==2, , drop=FALSE])  
  return(list(mcindex=mcindex))
}

#soft thresholding estimator of Huang et al. 2014+
# 3 lines added to handle flat case when largest total signal < p*bkgd
.sigclustcovest <- function(vsampeigv, sig2b) {
  d <- length(vsampeigv)
    #Check have some eigenvalues < sig2b
  vtaucand <- vsampeigv - sig2b
  #need to make sure there is enough total power.
  # if not, just return flat estimate as in Matlab impl, PKK 01/14/2014
  if (sum(vsampeigv) <= d*sig2b) {
    return(list(veigvest=rep(sig2b, d), tau=0))
  }
  #find threshold to preserve power
  which <- which(vtaucand<=0)
  icut <- which[1] - 1
  powertail <- sum(vsampeigv[(icut+1):d])
  power2shift <- sig2b*(d-icut) - powertail
  vi <- c(1:icut)
  vcumtaucand <- sort(cumsum(sort(vtaucand[vi])),decreasing=TRUE)
  vpowershifted <- (vi-1)*vtaucand[vi] + vcumtaucand
  flag <- vpowershifted < power2shift
  if(sum(flag)==0){
    itau <- 0;
  }else{
    which <- which(flag>0)
    itau <- which[1]
  }
  if(itau==1){
    powerprop <- power2shift/vpowershifted[1] #originally (../vpowershifted) no [1] idx
    tau <- powerprop*vtaucand[1]
  }else if(itau==0){
    powerprop <- power2shift/vpowershifted[icut] 
    tau <- powerprop*vtaucand[icut] 
  }else{
    powerprop <- (power2shift-vpowershifted[itau])/
      (vpowershifted[itau-1]-vpowershifted[itau]) 
    tau <- vtaucand[itau] + powerprop*(vtaucand[itau-1] - vtaucand[itau]) 
  }
  veigvest <- vsampeigv - tau 
  flag <- veigvest > sig2b 
  veigvest <- flag*veigvest + (1-flag)*(sig2b*rep(1,d))
  list(veigvest=veigvest, tau=tau)
}

#calculate variances of null Gaussian for simulation
# using factor model w/ eigenvalues of pca and background noise.
# copied from sigclust package
.vareigen <- function(x, n, p, icovest) {  
  #check the dimension of x to
  #match n and p
  if(dim(x)[1]==n & dim(x)[2]==p){
    mad1 <- mad(x)
    simbackvar <- mad1^2
    # Jan. 23, 07; replace eigen by
    #svd to save memory
    xcov <- cov(x)
    xeig <- eigen(xcov, symmetric=TRUE, only.values =TRUE)
    veigval <- xeig$values
    vsimeigval <- xeig$values
    
    #      avgx<-t(t(x)-colMeans(x))
    #	dv<-svd(avgx)$d
    #      veigval<-dv^2/(n-1)
    #	vsimeigval<-veigval
    
    if(icovest==1){
      vsimeigval <- .sigclustcovest(veigval, simbackvar)$veigvest
    }
    
    if(icovest==2){
      vsimeigval[veigval<0] <- 0
    }
    
    if(icovest==3){  # Use original background noise thresholded estimate
      # (from Liu, et al, JASA paper)
      vsimeigval[veigval<simbackvar] <- simbackvar
    }
    list(veigval=veigval, simbackvar=simbackvar, vsimeigval=vsimeigval)
  }else{
    print("Wrong size of matrix x!")
    return(0)
  }
}

#given null eigenvalues, simulate Gaussian datasets and compute
# cluster index - cleaned data simulation and changed to .simcluster(), PKK
.simnull <- function(vsimeigval, n, p, metric, linkage, square, l, nCIs) {
  simnorm <- matrix(rnorm(n*p, sd=sqrt(vsimeigval)), n, p, byrow=TRUE)
  simclust <- .simcluster(simnorm, p, metric, linkage, square, l, nCIs)
  list(mcindex=simclust$mcindex)
}







