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
#' @param K an integer value specifying the range of nodes for which to 
#'        create diagnostic plots (default = nrow(obj$p_norm)-1)
#' @param fname a character string specifying the name of the output file,
#'        see details for more information on default behavior depending on
#'        \code{length(K)} (default = NULL)
#' @param ci_idx a numeric value between 1 and \code{length(obj$in_args$ci)}
#'        specifying which cluster index to use for creating plots (default = 1)
#' @param pty a string specifying the plot type that should be created
#'        see details for more information on different plot types available
#'        (default = "all")
#' @param ... other parameters passed to plots
#' 
#' @return
#' Prints plots to \code{paste0(fname, ".pdf")} if \code{fname} is specified or \code{length(K)>1}.
#'
#' @details
#' If \code{K} is a single value, the default behavior is to output the diagnostic
#' plot to the current device. This behavior is overriden if \code{fname} is specified
#' by the user. If \code{K} is a set of values, than the diagnostic plots are
#' printed to \code{paste0(fname, ".pdf")}. If \code{fname} is not specified, a default value,
#' "shc_diagnostic", is used.
#'
#' The \code{pty} parameter accepts the following options:
#' \itemize{
#' \item{\code{"background"}}: {empirical distribution of marginal data used for background variance estimation}
#' \item{\code{"qq"}}: {qqplot for checking null Gaussian assumption}
#' \item{\code{"diag"}}: {screeplot for checking covariance factor model}
#' \item{\code{"pvalue"}}: {empirical distribution of simulated CIs}
#' \item{\code{"all"}}: {all of the above plots (default).}
#' }
#' 
#' 
#' @export
#' @name diagnostic-shc
#' @author Patrick Kimes
diagnostic.shc <- function(obj, K = nrow(obj$p_norm)-1,
                           fname = NULL, ci_idx = 1,
                           pty = "all", ...) {

    ## available plot types
    avail_pty <- c("background", "qq", "diag", "pvalue", "all")
    
    ## check default values
    if (length(K) == 0 || min(K) < 1 || max(K) > nrow(obj$p_norm)-1) {
        stop("K must be a set of indices between 1 and n-1")
    }
    if (!is.null(fname) || !is.character(fname)) {
        stop("fname must be NULL or a string")
    }
    if (length(pty) > 1) {
        stop("can only specify one plot type")
    } else if (sum(grepl(pty, avail_pty)) > 1) {
        stop("specified plot type matches more than one plot type")
    } else if (sum(grepl(pty, avail_pty)) == 0) {
        stop(paste("pty is not a valid value, please specify one of (or part of):",
                   paste(availty, collapse=", ")))
    }
    
    ## parse default fname behavior
    if (is.null(fname) && length(K) > 1) {
        fname <- "shc_diagnostic"
    }

    ## create plots 
    if (is.null(fname)) {
        .diagnostic_k(obj, K, ci_idx, pty)
    } else {
        pdf(paste0(fname, ".pdf"))
        for (k in K) {
            plot(0, 0, col="white")
            text(0, 0, paste("K =", k))
            .diagnostic_k(obj, k, ci_idx, pty)
        }
        dev.off()
    }
    
}



## code cleaned/modified from sigclust CRAN package, plot function
.diagnostic_k <- function(shc, k, ci_idx, pty) {

    ## identify subtree of x
    idx_sub <- unlist(shc$hc_idx[k, ])
    k_mat <- shc$in_mat[idx_sub, ]
    n <- nrow(k_mat)
    d <- ncol(k_mat)

    ## parameters from shc object
    icovest <- shc$in_args$icovest
    n_sim <- shc$in_args$n_sim
    keigval_dat <- shc$eigval_dat[k, ]
    backvar_k <- shc$backvar[k]
    keigval_sim <- shc$eigval_sim[k, ]
    kci_sim <- shc$ci_sim[k, , ci_idx]
    kci_dat <- shc$ci_dat[k, ci_idx]
    kp_emp <- shc$p_emp[k, ci_idx]
    kp_norm <- shc$p_norm[k, ci_idx]

    ## vectorize data and compute statistics
    overlay.x <- as.vector(k_mat)
    mean <- mean(overlay.x)
    sd <- sd(overlay.x)
    median <- median(overlay.x)
    mad <- mad(overlay.x)
    ntot <- length(overlay.x)
    maxnol <- 5000

    ## Background Standard Deviation Diagnostic Plot
    if (any(grepl(pty, c("all", "background")))) {
        par(mfrow=c(1, 1))
        par(mar=c(5, 4, 4, 2) + 0.1)
        nused <- maxnol
        denraw <- density(overlay.x)
        if (ntot > maxnol) {
            overlay.x <- overlay.x[sample(c(1:ntot), maxnol)]
        }else{
            nused <- ntot
        }
        
        xmin <- min(denraw$x)
        xmax <- max(denraw$x)
        ymin <- min(denraw$y)
        ymax <- max(denraw$y)
        
        overlay.y <- ymin + (0.15+0.5*runif(nused))*(ymax - ymin)
        plot(denraw, xlim=range(overlay.x), ylim=range(denraw$y),
             col="blue", xlab="", main="", lwd=3, ...)
        xgrid <- seq(xmin, xmax, by=0.0025*(xmax-xmin))
        normden <- dnorm(xgrid, mean=median, sd=sd)
        lines(xgrid, normden, col="red", lwd=3)
        points(overlay.x, overlay.y, col="green", pch=".")
        
        title(paste0("Distribution of All Pixel values combines, K=", k))
        if (ntot > maxnol) {
            text(xmin+0.47*(xmax-xmin), ymin+0.9*(ymax-ymin),
                 paste("Overlay of", as.character(maxnol), "of",
                       as.character(ntot), "data points"), cex=1.3)
        } else {
            text(xmin+0.47*(xmax-xmin), ymin+0.9*(ymax-ymin),
                 paste("Overlay of", as.character(ntot), "data points"),
                 cex=1.3)
        }
        text(xmin+0.47*(xmax-xmin), ymin+0.8*(ymax-ymin),
             paste("Mean =", as.character(round(mean, 3)),
                   "  Median =", as.character(round(median, 3))), cex=1.3)
        text(xmin+0.47*(xmax-xmin), ymin+0.7*(ymax-ymin),
             paste("s.d. =", as.character(round(sd, 3)),
                   "MAD =", as.character(round(mad, 3))), cex=1.3)
        text(xmin+0.47*(xmax-xmin), ymin+0.6*(ymax-ymin),
             paste0("Gaussian(", as.character(round(median, 3)), ",",
                    as.character(round(mad, 3)), ") density"),
             col="red", cex=1.3)
        if (mad > sd) {
            text(xmin+0.47*(xmax-xmin), ymin+0.55*(ymax-ymin),
                 "Warning: MAD > s.d., SHC can be anti-conservative",
                 cex=1.3)
        }
    }

    
    ## QQ plot
    if (any(grepl(pty, c("all", "qq")))) {
        par(mfrow=c(1, 1))  
        par(mar=c(5, 4, 4, 2) + 0.1)
        qqnorm <- qqnorm(as.vector(shc$in_mat), plot.it=FALSE)
        if(ntot > maxnol){
            which <- sample(c(1:ntot), maxnol)
        }else{
            which <- c(1:ntot)
        }
        x <- sort(qqnorm$x[which])
        y <- sort(qqnorm$y[which])
        x25 <- x[which(x > qnorm(0.25))[1]]
        x50 <- x[which(x > qnorm(0.5))[1]]
        x75 <- x[which(x > qnorm(0.75))[1]]
        y25 <- y[which(x > qnorm(0.25))[1]]
        y50 <- y[which(x > qnorm(0.5))[1]]
        y75 <- y[which(x > qnorm(0.75))[1]]
        plot(x, y, col="red", xlab="Gaussian Q", ylab="Data Q",
             main="Robust Fit Gaussian Q-Q, All Pixel values",
             cex.lab=1.3,...)
        abline(0, 1, col="green", lwd=2)
        xmin <- min(x) - 0.05*(max(x)-min(x))
        xmax <- max(x) + 0.05*(max(x)-min(x))
        ymin <- min(y) - 0.05*(max(y)-min(y))
        ymax <- max(y) + 0.05*(max(y)-min(y))
        text(xmin+0.3*(xmax-xmin), ymin+0.9*(ymax-ymin),
             paste("Mean =", as.character(round(mean, 3))), cex=1.3)
        text(xmin+0.3*(xmax-xmin), ymin+0.8*(ymax-ymin),
             paste("sd =", as.character(round(sd, 3))), cex=1.3)
        text(x25, y25, "+", cex=1.3)
        text(x25+0.7, y25, "0.25 quantile", cex=1.3)
        text(x50, y50, "+", cex=1.3)
        text(x50+0.7, y50, "0.5 quantile", cex=1.3)
        text(x75, y75, "+", cex=1.3)
        text(x75+0.7, y75, "0.75 quantile", cex=1.3)
    }


    ## covariance estimation diagnostic plot
    if (any(grepl(pty, c("all", "diag")))) {
        ncut <- 100
        if (d > ncut) {
            par(mfrow=c(2, 2))
        } else {
            par(mfrow=c(1, 2))
        }
        par(mar=c(2, 3.7, 2, 1.7))
        keigval_pos <- keigval_dat[which(keigval_dat > 10^(-12))]
        dpos <- length(keigval_pos)
        xmin <- 0
        xmax <- d+1
        ymin <- min(keigval_dat) - 0.05*(max(keigval_dat)-min(keigval_dat))
        ymax <- max(keigval_dat) + 0.05*(max(keigval_dat)-min(keigval_dat))
        plot(1:d, keigval_sim, type="l", lty=2, lwd=3, col="red",
             xlim=c(xmin, xmax), ylim=c(ymin, ymax),
             xlab="Component #", ylab="Eigenvalue", ...)
        points(1:d, keigval_dat, col="black")
        title(paste0("Eigenvalues, K=", k))
        lines(c(ncut+0.5, ncut+0.5), c(ymin, ymax), col="green")
        
        if (icovest != 2) {
            lines(c(0, d+1), c(backvar_k, backvar_k), col="magenta")
            text(xmin+0.45*(xmax-xmin), ymin+0.9*(ymax-ymin),
                 paste0("Background variance = ",
                        as.character(round(backvar_k, 3))),
                 col="magenta")
        }
        text(xmin+0.45*(xmax-xmin), ymin+0.8*(ymax-ymin),
             "Eigenvalues for simulation", col="red")
        if (mad > sd) {
            text(xmin + 0.45*(xmax-xmin), ymin + 0.65*(ymax-ymin),
                 "Warning: MAD > s.d.", col="magenta")
        }
        
        ymin <- min(log10(keigval_pos)) - 0.05*(max(log10(keigval_pos)) - min(log10(keigval_pos)))
        ymax <- max(log10(keigval_pos)) + 0.05*(max(log10(keigval_pos)) - min(log10(keigval_pos)))
        plot(1:d, log10(keigval_sim), type="l", lty=2, lwd=3, col="red",
             xlim=c(xmin, xmax), ylim=c(ymin, ymax),
             xlab="Component #", ylab="log10(Eigenvalue)", ...)
        points(1:dpos, log10(keigval_pos), col="black")
        title(paste0("log10 Eigenvalues, K=", k))
        lines(c(ncut+0.5, ncut+0.5), c(ymin, ymax), col="green")
        
        if (icovest != 2) {
            lines(c(0,d+1), log10(c(backvar_k, backvar_k)), col="magenta")
            text(xmin+0.45*(xmax-xmin), ymin+0.9*(ymax-ymin),
                 paste0("log10 Background variance = ",
                        as.character(round(log10(backvar_k), 3))),
                 col="magenta")
        }
        text(xmin+0.45*(xmax-xmin), ymin+0.8*(ymax-ymin),
             "Eigenvalues for simulation", col="red")
        if (mad > sd) {
            text(xmin+0.45*(xmax-xmin), ymin+0.65*(ymax-ymin),
                 "SHC may be Anti-iConservative", col="magenta")
        }
        if (length(keigval_dat) >= ncut) {
            xmin <- 0
            xmax <- ncut+1
            ymin <- min(keigval_dat[1:ncut]) - 0.05*(max(keigval_dat[1:ncut]) - min(keigval_dat[1:ncut])) 
            ymax <- max(keigval_dat[1:ncut]) + 0.05*(max(keigval_dat[1:ncut]) - min(keigval_dat[1:ncut])) 
            plot(1:ncut, keigval_sim[1:ncut], type="l", lty=2, lwd=3, col="red",
                 xlim=c(xmin, xmax), ylim=c(ymin, ymax),
                 xlab="Component #", ylab="Eigenvalue", ...)
            points(1:ncut, keigval_dat[1:ncut], col="black")
            title(paste0("Zoomed in version of above, K=", k))
            
            if (icovest != 2) {
                lines(c(0, d+1), c(backvar_k, backvar_k), col="magenta")
                text(xmin+0.45*(xmax-xmin), ymin+0.9*(ymax-ymin),
                     paste0("Background variance = ",
                            as.character(round(backvar_k, 3))),
                     col="magenta")
            }
            text(xmin+0.45*(xmax-xmin), ymin+0.8*(ymax-ymin),
                 "Eigenvalues for simulation", col="red")
            
            nmax <- min(dpos, ncut)
            ymin <- min(log10(keigval_pos[1:nmax])) - 0.05*(max(log10(keigval_pos[1:nmax])) - min(log10(keigval_pos[1:nmax])))
            ymax <- max(log10(keigval_pos[1:nmax])) + 0.05*(max(log10(keigval_pos[1:nmax])) - min(log10(keigval_pos[1:nmax])))
            plot(1:nmax,log10(keigval_sim[1:nmax]), type="l", lty=2, lwd=3,
                 col="red", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
                 xlab="Component #", ylab="log10(Eigenvalue)", ...)
            points(1:nmax, log10(keigval_pos[1:nmax]), col="black")
            title(paste0("log10 Eigenvalues, K=", k))
            
            if (icovest != 2) {
                lines(c(0, ncut+1), log10(c(backvar_k, backvar_k)), col="magenta")
                text(xmin+0.30*(xmax-xmin), ymin+0.9*(ymax-ymin),
                     paste0("log10 Background variance = ",
                            as.character(round(log10(backvar_k), 3))),
                     col="magenta")
            }
            text(xmin+0.30*(xmax-xmin), ymin+0.8*(ymax-ymin),
                 "Eigenvalues for simulation", col="red")
        }
    }

    
    ## p-value plot
    if (any(grepl(pty, c("all", "pvalue")))) {
        par(mfrow=c(1, 1))
        par(mar=c(5, 4, 4, 2) + 0.1)
        denpval <- density(kci_sim)
        denrange <- quantile(denpval$y, probs=c(0,0.5, 0.75, 1))
        dy <- 0.1*(denrange[4]-denrange[1])
        
        mindex <- mean(kci_sim)
        sindex <- sd(kci_sim)
        
        xmin <- min(c(kci_sim, kci_dat))
        xmax <- max(c(kci_sim, kci_dat))
        dx <- 0.1*(xmax-xmin)
        xind <- seq(xmin, xmax, 0.001)
        
        plot(denpval, xlim=c(xmin-dx, xmax+dx), col="red",
             xlab="Cluster Index", main="", lwd=2, ...)
        title(paste0("SHC Results, K=", k))
        points(kci_sim, runif(n_sim, denrange[2], denrange[3]),
               col="blue", pch=".", cex=2)
        lines(c(kci_dat, kci_dat), c(denrange[1]-dy, denrange[4])+dy,
              col="green", lty=2, lwd=2)
        lines(xind, dnorm(xind, mean=mindex, sd=sindex), col="black", lty=3, lwd=2)
        legend(xmin+0.05*(xmax-xmin), denrange[4], paste("P-value=", kp_emp),
               text.col="red", bty="n")
        legend(xmin+0.05*(xmax-xmin), denrange[3]+0.75*(denrange[4]-denrange[3]),
               paste("P-vNorm=", round(kp_norm, 3)), bty="n")
    }

}

