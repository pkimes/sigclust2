#' statistical Significance for Hierarchical Clustering (SHC) 
#' algorithm
#'
#' \code{shc} first computes the desired dendrogram, and
#' quantifies the statistical significance of clustering at each 
#' branching along the tree using a hypothesis testing procedure. 
#' The code is written so that various cluster indices can be 
#' easily introduced by writting new CI function and adding to 
#' .initcluster(), .simcluster(). When possible, this function 
#' makes use of a C++ implementation of hierarchical clustering 
#' available through the \code{Rclusterpp.hclust} package for the 
#' case of clustering by Pearson correlation (\code{dist="cor"}), 
#' we make use of \code{WGCNA::cor} with the usual 
#' code{stats::hclust}.
#' 
#' @param x a dataset with n rows and p columns, with observations in rows
#' @param metric a string specifying the metric to be used in the hierarchical 
#'        clustering procedure. This must be a metric accepted by \code{dist}, 
#'        e.g. "euclidean," or "cor."
#' @param linkage a string specifying the linkage to be used in the hierarchical 
#'        clustering procedure. This must be a linkage accepted by 
#'        \code{hclust}, e.g. "ward."
#' @param l an integer value specifying the power of the Minkowski distance, if 
#'        used, default is 2.
#' @param alphaStop a value between 0 and 1 specifying the desired level of the 
#'        test. If no FWER control is desired, simply set alpha to 1. The 
#'        default is 1. The testing procedure will terminate when no branches
#'        meet the corresponding FWER control threshold. This procedure may 
#'        substantially speed up the procedure by reducing the number of tests 
#'        considered. NOTE: if p-values are desired for all branches, the FWER
#'        cutoffs may be provided a posteriori by a call to \code{FWERcutoffs()}.
#'        Additionally, they may be specified at \code{plot()} using the 
#'        \code{alpha} and \code{FWER} parameters.
#' @param nsim a numeric value specifying the number of simulations for SigClust 
#'        testing. The default is to run 100 simulations at each merge. 
#' @param minObs an integer specifying the minimum number of observations needed
#'        to calculate a p-value. Default is 10.
#' @param icovest a numeric value specifying the covariance estimation method: 
#'        1. Use a soft threshold method as constrained MLE (default); 
#'        2. Use sample covariance estimate (recommended when diagnostics fail); 
#'        3. Use original background noise thresholded estimate (from Liu et 
#'        al., (2008)) ("hard thresholding") as described in the \code{sigclust}
#'        package documentation.
#' @param bkgdPCA a logical value whether to use principal component scores when
#'        estimating background noise under the null. Default is TRUE.
#' @param useCpp a logical value whether to use the \code{Rclusterpp} package.
#'        Default is TRUE.
#' @param nThreads a integer value specifying the number of threads to allocate 
#'        for the Rclusterpp process. Default is 1.
#' @param verb a logical value specifying whether the method should print out 
#'        when testing completes along each node along the dendrogram.
#'        Default is FALSE.
#' @param testCIs a string vector specifying the cluster indices to be used for 
#'        testing along the dendrogram. Currently, options include: "2CI", 
#'        "linkage". Default is "2CI". 
#' @param testNulls a string vector specifying the clustering approach that 
#'        should be used at each node as the comparison. Currently, options
#'        include: "2means", "hclust". Note, testNulls and testCIs must be of
#'        equal length. Default is "hclust".
#' @param cutoffCI a single value between 1 and \code{length(testCIs)} 
#'        specifiying which CI to use for the FWER stopping rule.
#'        This only has an effect if alphaStop is specified to a non-default 
#'        value. Default is 1.
#' 
#' @return The function returns a \code{shc} object containing the 
#'         resulting p-values. The print method call will output a dendrogram
#'         with the corresponding p-values placed at each merge. 
#' 
#' @details The function expands on the \code{sigclust} idea to the hierarchical 
#'          setting by modifying the clustering procedure employed to compute 
#'          the null distribution of cluster indices.
#' 
#' @examples
#' hsc_cars <- shc(mtcars, metric="euclidean", linkage="single")
#' tail(mpvalnorm(hsc_cars), 10)
#' 
#' @import Rclusterpp WGCNA
#' @export shc
#' @author Patrick Kimes


shc <- function(x, metric, linkage, alphaStop=1, l=2, bkgdPCA=TRUE,
                nsim=100, minObs=10, icovest=1, useCpp=TRUE, nThreads=1,
                testCIs="2CI", testNulls="hclust", cutoffCI=1) {  

    if (useCpp) {
        Rclusterpp.setThreads(nThreads)
    }
    
    ##number of cluster indices
    nCIs <- length(testCIs)
    ##check validity of testCIs/testNulls
    if (length(testCIs) != length(testNulls)) {
        nCIs <- min(length(testCIs), length(testNulls))
        cat("!! testCIs and testNulls must be of equal length.  !!")
        cat(paste("!! Only using the first", nCIs, "entries of each.         !!"))
    }
    
    ##check validity of cutoffCI
    if (cutoffCI > nCIs) { cutoffCI <- 1 }
    
    ##check validity of alphaStop
    if (alphaStop > 1 || alphaStop < 0) { alphaStop <- 1}
    
    
    ##check the dimension of x to match n and p
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    x <- as.matrix(x)
    xclust <- .initcluster(x, n, p, metric, linkage, l, 
                           nCIs, testCIs, useCpp)
    xmcindex <- xclust$mcindex
    
    ##need to correct height of merges if using Ward clustering w/ Rclusterpp
    if (linkage=="ward") {
        xclust$clusters$height <- sqrt(2*xclust$clusters$height)
    }
    
    ##p-values for all <=(n-1) tests
    mpval <- matrix(2, nrow=n-1, ncol=nCIs)
    mpvalnorm <- matrix(2, nrow=n-1, ncol=nCIs)
    colnames(mpval) <- paste(testNulls, testCIs, sep="_")
    colnames(mpvalnorm) <- paste(testNulls, testCIs, sep="_")
    
    ##null covariance parameters for all <=(n-1) tests
    meigval <- matrix(-1, nrow=n-1, ncol=p)
    msimeigval <- matrix(-1, nrow=n-1, ncol=p)
    vsimbackvar <- rep(-1, n-1)
    asimcindex <- array(-1, dim=c(n-1, nsim, nCIs))

    ##determine parent branch node for all children nodes along dendrogram
    allPDpairs <- rbind(cbind(xclust$clusters$merge[,1], 1:(n-1)), 
                        cbind(xclust$clusters$merge[,2], 1:(n-1)))
    PDmap <- data.frame(allPDpairs[allPDpairs[, 1]>0, ])
    names(PDmap) <- c("daughter", "parent")
    PDmap <- PDmap[order(PDmap$daughter), 2] #the parent of each daughter
    PDmap <- c(PDmap, n) #add final node without a parent
    
    ##compute Meinshausen cutoffs for significance at alpha
    clusterSizes <- apply(xclust$clusterList, 1, 
                          function(x) { length(unlist(x)) })
    cutoff <- alphaStop * clusterSizes/n
    
    treesig <- rep(TRUE, n)
    
    for (k in (n-1):1) { #move through layers in reverse order.
        
        ##if parent wasn't significant, skip
        if ( !treesig[PDmap[k]] ) {
            mpvalnorm[k, ] <- rep(47, nCIs)
            mpval[k, ] <- rep(47, nCIs)
            meigval[k, ] <- rep(-47, p)
            msimeigval[k, ] <- rep(-47, p)
            vsimbackvar[k] <- -47
            
            treesig[k] <- FALSE
            next 
        }
        
        ## THESE CALCULATIONS CAN PROBABLY BE ROLLED INTO
        ## CLUSTERSIZES/CUTOFF COMPUTATION
        subxIdx <- unlist(xclust$clusterList[k, ])
        subn <- length(subxIdx)
        
        ##only calc p-values for branches w/ more than minObs
        if (subn >= minObs) {
            xvareigen <- .vareigen(x[subxIdx, ], subn, p, icovest, bkgdPCA)
            for (i in 1:nsim) {
                xsim <- .simnull(xvareigen$vsimeigval, subn, p, 
                                 metric, linkage, l, 
                                 nCIs, testCIs, testNulls, useCpp)
                asimcindex[k, i, ] <- xsim$mcindex        
            }
            mindex <- colMeans(as.matrix(asimcindex[k, , ]))
            sindex <- apply(as.matrix(asimcindex[k, , ]), 2, sd)
            mpvalnorm[k, ] <- pnorm(xmcindex[k, ], mindex, sindex)
            mpvalnorm[k, testCIs=="linkage"] <- 
                1-mpvalnorm[k, testCIs=="linkage"] #flip for linkage based testing
            mpval[k, ] <- colMeans(as.matrix(asimcindex[k, , ]) <= 
                                       matrix(xmcindex[k, ], nrow=nsim, 
                                              ncol=nCIs, byrow=TRUE))
            mpval[k, testCIs=="linkage"] <-
                1-mpval[k, testCIs=="linkage"] #flip for linkage based testing
            meigval[k, ] <- xvareigen$veigval
            msimeigval[k, ] <- xvareigen$vsimeigval
            vsimbackvar[k] <- xvareigen$simbackvar
            
            ##update treesig based on significance threshold if 
            ## short-stopping rule is desired
            if (alphaStop < 1) {
                treesig[k] <- (mpvalnorm[k, cutoffCI] < cutoff[k])        
            }
            
        }    
        
    }  
    
    return(new("shc", 
               data = x,
               meigval = meigval,
               msimeigval = msimeigval,
               vsimbackvar = vsimbackvar,
               asimcindex = asimcindex,
               mpval = mpval,
               mpvalnorm = mpvalnorm,
               xmcindex = xmcindex,
               clusters = xclust$clusterList,
               hc = xclust$clusters,
               args = list(metric = metric,
                   linkage = linkage,
                   alphaStop = alphaStop,
                   l = l,
                   bkgdPCA = bkgdPCA,
                   nsim = nsim,
                   minObs = minObs,
                   icovest = icovest,
                   testCIs = testCIs,
                   testNulls = testNulls,
                   cutoffCI = cutoffCI)))
}



##calculate 2-means cluster index (nxp matrices)
.sumsq <- function(x) { norm(sweep(x, 2, colMeans(x), "-"), "F")^2 }
.calc2CI <- function(x1, x2) {
    if (is.matrix(x1) && is.matrix(x2) && ncol(x1) == ncol(x2)) {
        (.sumsq(x1) + .sumsq(x2)) / .sumsq(rbind(x1, x2))
    } else {
        error(paste("x1, x2 must be matrices with same ncols",
                    "for 2CI calculation"))
    }      
}



##perform hierarchical clustering on the original data and 
## compute the corresponding cluster indices for each merge
.initcluster <- function(x, n, p, metric, linkage, l, 
                         nCIs, testCIs, useCpp) { 

    ##need to implement clustering algorithm
    if (metric == 'cor') {
        dmatrix <- 1 - WGCNA::cor(t(x))
        clusters <- hclust(as.dist(dmatrix), method=linkage)
    } else {
        if (useCpp) {
            clusters <- Rclusterpp.hclust(x, method=linkage, 
                                          distance=metric, p=l)
        } else {
            clusters <- hclust(dist(x, method=metric, p=l), method=linkage)
        }
    }

    ##matrix containing cluster indices
    mcindex <- matrix(-1, nrow=n-1, ncol=nCIs)
    
    ##list array of cluster indices at each of the n-1 merges
    clusterList <- array(list(), c(2*n-1, 2))
    clusterList[1:n, 1] <- as.list(n:1)
    clusterList[(n+1):(2*n-1), ] <- clusters$merge+n+(clusters$merge<0)

    ##complete clusterList
    for (k in 1:(n-1)) {
        clusterList[[n+k, 1]] <- unlist(clusterList[clusterList[[n+k, 1]], ])
        clusterList[[n+k, 2]] <- unlist(clusterList[clusterList[[n+k, 2]], ])
    }
    
    ##calculate cluster index(ices) for merge k
    for (iCI in 1:nCIs) {
        if (testCIs[iCI] == "2CI") {
            for (k in 1:(n-1)) {
                mcindex[k, iCI] <- .calc2CI(x[clusterList[[n+k, 1]], , drop=FALSE],
                                            x[clusterList[[n+k, 2]], , drop=FALSE])
            }
        } else if (testCIs[iCI] == "linkage") {
            mcindex[, iCI] <- clusters$height
        }
    }
    clusterList <- clusterList[-(1:n), ]
    
    list(clusters=clusters, 
         clusterList=clusterList,
         mcindex=mcindex))
}



##given null eigenvalues, simulate Gaussian datasets and compute
## cluster index - cleaned data simulation and changed to .simcluster(), PKK
.simnull <- function(vsimeigval, n, p, metric, linkage, l, 
                     nCIs, testCIs, testNulls, useCpp) {
    simnorm <- matrix(rnorm(n*p, sd=sqrt(vsimeigval)), 
                      n, p, byrow=TRUE)
    simclust <- .simcluster(simnorm, p, metric, linkage, l, 
                            nCIs, testCIs, testNulls, useCpp)
    list(mcindex=simclust$mcindex)
}



##perform hierarchical clustering on a simulated dataset and
## compute the correspond cluster indices for only the final merge
.simcluster <- function(sim_x, p, metric, linkage, l, 
                        nCIs, testCIs, testNulls, useCpp) { 
    ##need to implement clustering algorithm
    if (metric == "cor") {
        dmatrix <- 1 - WGCNA::cor(t(sim_x))
        clusters <- hclust(as.dist(dmatrix), method=linkage)
    } else {
        if (useCpp) {
            clusters <- Rclusterpp.hclust(sim_x, method=linkage, 
                                          distance=metric, p=l)      
        } else {
            clusters <- hclust(dist(sim_x, method=metric, p=l), method=linkage)      
        }
    }
    sim_split <- cutree(clusters, k=2)
    mcindex <- matrix(-1, nrow=1, ncol=nCIs)
    for (iCI in 1:nCIs) {
        if (testCIs[iCI] == "2CI") {
            if (testNulls[iCI] == "hclust") {
                mcindex[iCI] <- .calc2CI(sim_x[sim_split==1, , drop=FALSE],
                                         sim_x[sim_split==2, , drop=FALSE])        
            } else if (testNulls[iCI] == "2means") {
                kmsol <- kmeans(sim_x, centers=2)
                mcindex[iCI] <- kmsol$tot.withinss/kmsol$totss
            }
        } else if (testCIs[iCI] == "linkage") {
            mcindex[iCI] <- clusters$height[nrow(sim_x)-1]
        }
    }

    list(mcindex=mcindex)
}









