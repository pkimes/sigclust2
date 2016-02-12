#' plot diagnostics for shc object
#'
#' Provides visualizations to check the null Gaussian assumption at a specified
#' range of nodes along the dendrogram
#' 
#' @param obj a \code{shc} object
#' @param K an integer value or vector specifying the range of nodes for which to 
#'        create diagnostic plots, with 1 corresponding to the root (default = 1)
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
#' by the user. If \code{length(K) > 1}, than the diagnostic plots are
#' printed to \code{paste0(fname, ".pdf")}. If \code{fname} is not specified,
#' "shc_diagnostic", is used by default.
#'
#' The \code{pty} parameter accepts the following options:
#' \itemize{
#' \item{\code{"background"}}: {empirical distribution of marginal data used for background variance estimation}
#' \item{\code{"qq"}}: {qqplot for checking null Gaussian assumption}
#' \item{\code{"covest"}}: {screeplot for checking covariance factor model}
#' \item{\code{"pvalue"}}: {empirical distribution of simulated CIs}
#' \item{\code{"all"}}: {all of the above plots (default).}
#' }
#' 
#' @import ggplot2 ggthemes
#' @name diagnostic-shc
#' @export
#' @method diagnostic shc
#' @author Patrick Kimes
diagnostic.shc <- function(obj, K = 1, fname = NULL, ci_idx = 1,
                           pty = "all", ...) {

    ## available plot types
    avail_pty <- c("background", "qq", "covest", "pvalue", "all")
    
    ## check default values
    if (length(K) == 0 || min(K) < 1 || max(K) > nrow(obj$p_norm)-1) {
        stop("K must be a set of indices between 1 and n-1")
    }
    if (!is.null(fname) && !is.character(fname)) {
        stop("fname must be NULL or a string")
    }
    if (length(pty) > 1) {
        stop("can only specify one plot type")
    } else if (sum(grepl(pty, avail_pty)) > 1) {
        stop("specified plot type matches more than one plot type")
    } else if (sum(grepl(pty, avail_pty)) == 0) {
        stop(paste("pty is not a valid value, please specify one of (or part of):",
                   paste(avail_pty, collapse=", ")))
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
.diagnostic_k <- function(obj, k, ci_idx, pty) {

    ## identify subtree of x
    idx_sub <- unlist(obj$idx_hc[k, ])
    k_mat <- obj$in_mat[idx_sub, ]
    n <- nrow(k_mat)
    d <- ncol(k_mat)

    ## parameters from obj object
    icovest <- obj$in_args$icovest
    n_sim <- obj$in_args$n_sim
    keigval_dat <- obj$eigval_dat[k, ]
    backvar_k <- obj$backvar[k]
    keigval_sim <- obj$eigval_sim[k, ]
    kci_sim <- obj$ci_sim[k, , ci_idx]
    kci_dat <- obj$ci_dat[k, ci_idx]
    kp_emp <- obj$p_emp[k, ci_idx]
    kp_norm <- obj$p_norm[k, ci_idx]

    ## vectorize data and compute statistics
    xvec <- as.vector(k_mat)
    mean_x <- mean(xvec)
    sd_x <- sd(xvec)
    med_x <- median(xvec)
    mad_x <- mad(xvec)
    ntot <- length(xvec)
    max_nol <- 5000


    ## Background Standard Deviation Diagnostic Plot
    if (any(grepl(pty, c("all", "background")))) {

        ## fit density to full data
        den_x <- density(xvec)
        den_df <- data.frame(x=den_x$x, y=den_x$y)

        ## subset on dataset for plotting
        nused <- min(max_nol, ntot)
        if (ntot > max_nol) {
            xvec <- xvec[sample(1:ntot, max_nol)]
        }

        ## determine range of density plot
        xmin <- min(den_df$x)
        xmax <- max(den_df$x)
        ymin <- min(den_df$y)
        ymax <- max(den_df$y)

        ## create df with plotting points
        xvec_df <- data.frame(x=xvec, y=runif(length(xvec), ymax*1/4, ymax*3/4))

        ## plot data kde and best fit gaussian
        gp <- ggplot(xvec_df) +
            geom_point(aes(x=x, y=y), alpha=1/5, size=1, color="#33a02c") +
            geom_path(aes(x=x, y=y), data=den_df, color="black", size=1) + 
            stat_function(fun=dnorm, args=list(mean=med_x, sd=mad_x),
                          color="#1f78b4", size=1, n=500) +
            theme_gdocs() +
            ylab("density") +
            ggtitle(paste0("Distribution of Vectorized Data for K=", k)) +
                scale_x_continuous(expand=c(0.01, 0)) +
                    scale_y_continuous(expand=c(0.01, 0))

        ## note # sample points
        if (ntot > max_nol) {
            gp <- gp +
                annotate(geom="text", x=xmax, y=ymax*.98, vjust=1, hjust=1,
                         label=paste("overlay of", max_nol, "/", ntot, "data points"),
                         color="#33a02c")
        } else {
            gp <- gp +
                annotate(geom="text", x=xmax, y=ymax*.98, vjust=1, hjust=1,
                         label=paste("overlay of", ntot, "data points"),
                         color="#33a02c")
        }

        ## label best fit gaussian density
        gp <- gp +
            annotate(geom="text", x=xmax, y=ymax*.90, vjust=1, hjust=1,
                     label=paste0("N(", round(med_x, 3), ", ", round(mad_x, 3), ") density"),
                     color="#1f78b4")

        ## return parameter values
        gp <- gp +
            annotate(geom="text", x=xmin, y=ymax*.98, vjust=1, hjust=0,
                     label=paste0("mean = ", round(mean_x, 3), "\n",
                         "median = ", round(med_x, 3), "\n",
                         "s.d. = ", round(sd_x, 3), "\n",
                         "MAD = ", round(mad_x, 3)))

        ## check for possible anti-conservative behavior
        if ((mad_x > sd_x) && (icovest == 1)) {
            gp <- gp +
                annotate("text", x=(xmax+xmin)/2, y=ymax/4, vjust=1,
                         label="Warning: MAD > s.d., SHC can be anti-conservative",
                         size=5, color="#e41a1c")
        }
        print(gp)
    }
    
    
    ## QQ plot
    if (any(grepl(pty, c("all", "qq")))) {

        ## compute quantiles
        qq_x <- qqnorm(as.vector(obj$in_mat), plot.it=FALSE)
        qq_df <- data.frame(x=qq_x$x, y=qq_x$y)

        ## subset for plotting
        if(ntot > max_nol) {
            qq_df <- qq_df[sort(sample(1:ntot, max_nol)), ]
        }

        ## build base of plot
        gp <- ggplot(qq_df) +
            geom_abline(aes(intercept=0, slope=1), color="#33a02c", size=1/2) +
        geom_point(aes(x=x, y=y), size=1, color="#e31a1c", alpha=1/2) +
        ggtitle(paste0("Robust Fit Gaussian Q-Q for K=", k)) +
            xlab("Gaussian Q") + ylab("Data Q") +
                theme_gdocs()

        ## include annotations in plot
        gp <- gp +
            annotate(geom="text", vjust=1, hjust=0,
                     x=min(qq_df$x), y=max(qq_df$y),
                     label=paste0("mean = ", round(mean_x, 3), "\n",
                         "s.d. =", round(sd_x, 3)))
        gp <- gp + annotate(geom="point", x=qnorm(c(.25, .50, .75)),
                            y=quantile(qq_df$y, c(.25, .50, .75)),
                            color="#1f78b4", shape=43, size=8)
        gp <- gp + annotate(geom="text", vjust=1, hjust=0, x=qnorm(c(.25, .50, .75)),
                            y=quantile(qq_df$y, c(.25, .50, .75)),
                            label=paste(" ", c("0.25", "0.50", "0.75"), "quantile"))
        
        print(gp)
    }
    

    
    ## covariance estimation diagnostic plot
    if (any(grepl(pty, c("all", "covest")))) {
        
        ymin <- min(keigval_dat) - 0.05*(max(keigval_dat)-min(keigval_dat))
        ymax <- max(keigval_dat) + 0.05*(max(keigval_dat)-min(keigval_dat))

        df <- data.frame(idx=rep(1:d, times=2),
                         grp=rep(c("sim", "dat"), each=d),
                         eigval=c(keigval_sim, keigval_dat))

        ## create base plot
        gp <- ggplot(df) +
            geom_point(aes(x=idx, y=eigval, color=grp)) +
                geom_path(aes(x=idx, y=eigval, group=grp), alpha=1/2) +
                    xlab("Component #") + ylab("Eigenvalue") +
        ggtitle(paste0("Eigenvalues, K=", k)) +
        theme_gdocs() +
            scale_color_manual(limits=c("sim", "dat"),
                               values=c("#e41a1c", "black"), guide=FALSE)

        if (icovest != 2) {
            ## include horizontal line at bkgd noise level
            gp <- gp +
                geom_hline(yintercept=backvar_k, col="#377eb8") +
                    annotate("text", hjust=1, vjust=0,
                             x=d+1, y=backvar_k+(ymax-ymin)*.02,
                             label=paste0("bkgd var = ",
                                 round(backvar_k, 3)),
                             col="#377eb8")
        }

        ## annotate eigvals
        gp <- gp +
            annotate("text", hjust=1, vjust=1,
                     x=d+1, y=ymin+0.9*(ymax-ymin),
                     label="eigenvalues for simulation", col="#e41a1c") +
        annotate("text", hjust=1, vjust=1,
                 x=d+1, y=ymin+0.84*(ymax-ymin),
                 label="sample eigenvalues")

        ## make note if anti-conservative behavior expected 
        if (mad_x > sd_x) {
            gp <- gp + 
                annotate("text", hjust=1, vjust=1,
                         x=d+1, y=ymin + 0.65*(ymax-ymin),
                         label="Warning: MAD > s.d.", col="#377eb8")
        }
        print(gp)
    }
    

    
    ## p-value plot
    if (any(grepl(pty, c("all", "pvalue")))) {

        ## fit density to cluster indices
        den_pv <- density(kci_sim)
        den_df <- data.frame(x=den_pv$x, y=den_pv$y)
        
        ## determine range of density plot
        xmin <- min(den_df$x)
        ymax <- max(den_df$y)

        mindex <- mean(kci_sim)
        sindex <- sd(kci_sim)

        ## create df with plotting points
        kci_df <- data.frame(x=kci_sim, y=runif(length(kci_sim), ymax*1/4, ymax*3/4))

        ## plot data kde and best fit gaussian
        gp <- ggplot(kci_df) +
            geom_point(aes(x=x, y=y), alpha=1/2, size=1, color="#377eb8") +
            geom_path(aes(x=x, y=y), data=den_df, color="#e41a1c", size=1) + 
            stat_function(fun=dnorm, args=list(mean=mindex, sd=sindex),
                          color="black", size=1, n=500) +
            geom_vline(xintercept=kci_dat, color="#4daf4a", size=1) + 
            theme_gdocs() +
            ylab("density") + xlab("Cluster Index") + 
                ggtitle(paste0("SHC Results for ci_idx=", ci_idx, ", K=", k)) +
                    scale_x_continuous(expand=c(0.01, 0)) +
                        scale_y_continuous(expand=c(0.01, 0))

        ## return p-values
        gp <- gp +
            annotate(geom="text", x=min(xmin, kci_dat), y=ymax*.98, vjust=1, hjust=0,
                     label=paste0(" p-value (Q) = ", round(kp_emp, 3), "\n",
                         " p-value (Z) = ", round(kp_norm, 3)))
        print(gp)
    }

}

