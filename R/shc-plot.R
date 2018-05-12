#' plot shc object
#'
#' Visualize the results of SHC analysis as an annotated 
#' dendrogram with significant branches highlighted
#' 
#' @param x a \code{shc} object to plot produced by a call to \code{shc}
#' @param groups a vector specifying group labels for the clustered objects.
#'        The vector should be in the same order as the rows of the original 
#'        data matrix. If specified, color blocks will be placed along the 
#'        bottom of the dendrogram. Useful when the samples have a priori known
#'        groupi behavior. (default = \code{NULL})
#' @param use_labs a boolean specifyin whether rowlabels should be added as 
#'        text along the bottom of the dendrogram (default = \code{TRUE})
#' @param fwer a boolean specifying whether the FWER control procedure of 
#'        Meinshausen et al. 2010 should be used. (default = \code{TRUE})
#' @param alpha a double between 0 and 1 specifying the significance cutoff. If 
#'        \code{fwer} is TRUE, the FWER of the entire dendrogram is controlled at 
#'        \code{alpha}, else, each branch is tested at \code{alpha}. Only has
#'        effect if \code{alpha(shc)} = 1. (default = 0.05 or
#'        \code{alpha} used for \code{x} if original \code{shc} called
#'        with \code{alpha} < 1)
#' @param ci_idx a numeric value between 1 and \code{length(ci)} 
#'        specifiying which CI to use for the FWER stopping rule.
#'        This only has an effect if \code{alpha} < 1. (default = 1 or
#'        \code{ci_idx} used for \code{x} if original \code{shc} called with
#'        \code{alpha} < 1 or \code{fwer} = TRUE)
#' @param ci_emp a logical value specifying whether to use the empirical
#'        p-value from the CI based on \code{ci_idx} for the FWER stopping rule.
#'        As with \code{ci_idx} this only has an effect if \code{alpha} < 1.
#'        (default = FALSE or \code{ci_emp} used for \code{x} if original \code{shc}
#'        called with \code{alpha} < 1 or \code{fwer} = TRUE)
#' @param hang a double value corresponding to the \code{hang} parameter for 
#'        the typical call to \code{plot} for an object of class
#'        \code{hsigclust} (default = -1)
#' @param ... other parameters to be used by the function
#'
#' @return
#' \code{ggplot} object containing a dendrogram annotated by the results of the
#' corresponding \code{shc} analysis
#' 
#' @details
#' This function makes use of dendrogram plotting functions made
#' available through the \pkg{ggdendro} package which provides a
#' \pkg{ggplot2}-like grammer for working with dendrograms.
#' 
#' @import ggplot2 ggdendro dplyr
#' @importFrom stats as.dendrogram hclust is.leaf
#' @name plot-shc
#' @export
#' @method plot shc
#' @author Patrick Kimes
plot.shc <- function(x, groups = NULL, use_labs = TRUE, hang = -1,
                     fwer = TRUE, alpha = NULL, ci_idx = NULL, ci_emp = NULL,
                     ...) {
    shc <- x

    ## determine number of samples
    n <- nrow(shc$p_emp) + 1
    
    ## colors to be used in figure
    col_tx_sig <- "#A60A3C"
    col_nd_null <- "gray"
    
    ## check validity of ci_idx
    if (!is.null(ci_idx)) {
        if (ci_idx > ncol(shc$p_emp)) {
            stop("invalid choice for ci_idx; ci_idx must be < length(ci)")
        }
    }
    
    ## check validity of alpha
    if (!is.null(alpha)) {
        if (alpha > 1 || alpha < 0) {
            stop("invalid choice for alpha; alpha must be 0 < alpha < 1")
        }
    }

    ## check fwer, alpha, ci_idx, ci_emp feasibility w/ original alpha value
    if (shc$in_args$alpha < 1) {
        if (is.null(alpha)) {
            alpha <- shc$in_args$alpha
        } else if (alpha > shc$in_args$alpha) {
            alpha <- shc$in_args$alpha
            warning(paste0("shc constructed using smaller alpha, using alpha = ", alpha))
        }
        
        if (is.null(ci_emp)) {
            ci_emp <- shc$in_args$ci_emp
        } else if (shc$in_args$ci_emp != ci_emp) {
            ci_emp <- shc$in_args$ci_emp
            warning("shc constructed using alpha < 1, using ci_emp from x$in_args")
        }
        
        if (is.null(ci_idx)) {
            ci_idx <- shc$in_args$ci_idx
        } else if (shc$in_args$ci_idx != ci_idx) {
            ci_idx <- shc$in_args$ci_idx
            warning("shc constructed using alpha < 1, using ci_idx from x$in_args")
        }

        if (!fwer) {
            fwer <- TRUE
            warning("shc constructed using alpha < 1, using fwer = TRUE")
        }
        
    } else {
        ## set default values if shc called with alpha = 1
        if (is.null(alpha)) {
            alpha <- 0.05
        }
        if (is.null(ci_idx)) {
            ci_idx <- 1
        }
        if (is.null(ci_emp)) {
            ci_emp <- FALSE
        }
    }
    
    ## for easier calling
    if (ci_emp) {
        p_use <- shc$p_emp[, ci_idx]
    } else {
        p_use <- shc$p_norm[, ci_idx]
    }
    
    ## determine sig/non-sig nodes based on cutoffs
    if (fwer) {
        ## use FWER controlling hierarchical testing procedure
        cutoff <- fwer_cutoff(shc, alpha)
        pd_map <- .pd_map(shc$hc_dat, n)
        nd_type <- rep("", n-1)
        for (k in 1:(n-1)) {
            ## check if subtree is large enough
            if (length(unlist(shc$idx_hc[k, ])) < shc$in_args$n_min) {
                nd_type[k] <- "n_small"
                next
            }
            ## check if parent was significant
            if ((k > 1) && (nd_type[pd_map[k]] != "sig")) {
                nd_type[k] <- "no_test"
                next
            }
            ## compare against p-value cutoff
            if (alpha < 1) {
                nd_type[k] <- ifelse(p_use[k] < cutoff[k], "sig", "not_sig")
            }
        }

    } else {
        ##if not using FWER control, use alpha as flat cutoff
        cutoff <- rep(alpha, n-1)
        nd_type <- ifelse(p_use < alpha, "sig", "not_sig")
        nd_type[shc$nd_type == "n_small"] <- "n_small"
    }
        

    ##using ggdendro package
    shc_dend <- as.dendrogram(shc$hc_dat, hang=hang)
    shc_dendat <- ggdendro::dendro_data(shc_dend)
    shc_segs <- ggdendro::segment(shc_dendat)
    shc_labs <- ggdendro::label(shc_dendat)

    ##change label hight to be correct
    shc_labs$y <- .lab_height(shc_dend)
    
    if (!is.null(groups) & length(groups) == n) {
        shc_labs$clusters <- groups[shc$hc_dat$order]
    }

    ## flip index, hclust/dend and shc use revered ordering
    rev_nd_type <- rev(nd_type)
    rev_p_use <- rev(p_use)
    
    ## significant nodes
    sig_linkvals <- as.factor(shc$hc_dat$height[rev_nd_type == "sig"])
    sig_segs <- filter(shc_segs, as.factor(y) %in% sig_linkvals)
    
    ## non-signficant nodes
    notsig_linkvals <- as.factor(shc$hc_dat$height[rev_nd_type == "not_sig"])
    notsig_segs <- filter(shc_segs, as.factor(y) %in% notsig_linkvals)

    ## tested nodes
    test_linkvals <- as.factor(shc$hc_dat$height[which(rev_nd_type == "sig")])
    test_segs <- filter(shc_segs, as.factor(y) %in% test_linkvals)
    test_segtops <- filter(test_segs, (y == yend) & (x < xend))
    test_segtops <- test_segtops[order(test_segtops[, 2]), ]
    test_segtops <- cbind(test_segtops, 
                          "pval"=format(rev_p_use[which(rev_nd_type == "sig")], 
                                        digits=3, scientific=TRUE))
    
    ## small sample nodes
    skip_linkvals <- as.factor(shc$hc_dat$height[rev_nd_type == "n_small"])
    skip_segs <- filter(shc_segs, as.factor(y) %in% skip_linkvals)
    
    ## un-tested nodes (non-significant parent node)
    if (fwer) {
        fwer_linkvals <- as.factor(shc$hc_dat$height[rev_nd_type == "no_test"])
        fwer_segs <- filter(shc_segs, as.factor(y) %in% fwer_linkvals)
    }
        
    ## calculate various plotting dimensions prior to actual ggplot call
    ax_x_ref <- max(shc$hc_dat$height)
    ax_x_top <- max(ax_x_ref)*1.25
    ax_x_bot <- -max(ax_x_ref)/4
    ax_x_scale <- floor(log10(ax_x_top))
    ax_y_range <- max(shc_segs$y) - min(shc_segs$y)  
    
    
    ## make initial ggdendro with null color
    plot_dend <- ggplot() + 
        geom_segment(data=shc_segs, 
                     aes(x=x, y=y, xend=xend, yend=yend), 
                     color=col_nd_null) +
        theme_bw()
    
    ## add appropriate title
    plot_dend <- plot_dend + 
        ggtitle(paste("showing all p-values below", 
                      alpha, "cutoff", 
                      ifelse(fwer, "(FWER corrected)", "")))
    
    ##add color group labels for objects if specified
    if (!is.null(groups)) {
        shc_labs$color_y <- shc_labs$y - ax_y_range/30
        plot_dend <- plot_dend +
            geom_tile(data=shc_labs,
                      aes(x=x, y=color_y, fill=as.factor(clusters)),
                      height=ax_y_range/30, alpha=1) +
            scale_fill_discrete('Labels')
    }
    
    ##add text labels for each object if desired
    if (use_labs) {
        shc_labs$txt_y <- shc_labs$y - (2-is.null(groups))*ax_y_range/30
        plot_dend <- plot_dend +
            geom_text(data=shc_labs, 
                      aes(x=x, y=txt_y, label=label, 
                          hjust=1, vjust=.5, angle=90), 
                      size=3)
    }
    
    ##add colored segments to plot if any branches were significant
    if (sum(rev_nd_type == "sig") > 0) {
        plot_dend <- plot_dend +
            geom_segment(data=sig_segs,
                         aes(x=x, y=y, xend=xend, yend=yend, color="sig"), 
                         size=1)
    }

    
    ##add p-values for segments if they were tested
    if (sum(rev_nd_type == "sig") > 0) {
        plot_dend <- plot_dend +
            geom_text(data=test_segtops, 
                      aes(x=x, y=y, label=pval, hjust=-0.2, vjust=-0.5),
                      col=col_tx_sig, size=4)
    }
    
    
    ##add FWER controlled segments if any branches were controlled
    if (fwer & sum(rev_nd_type == "no_test") > 0) {
        plot_dend <- plot_dend +
            geom_segment(data=fwer_segs,
                         aes(x=x, y=y, xend=xend, yend=yend, color="no_test"), 
                         size=1)
    }

    ##add colored segments if any branches had size < n_min
    if (sum(rev_nd_type == "n_small") > 0) {
        plot_dend <- plot_dend +
            geom_segment(data=skip_segs,
                         aes(x=x, y=y, xend=xend, yend=yend, color="n_small"), 
                         size=1)
    }
    
    ## add colored segments if any branches had non-significant calls
    if (sum(rev_nd_type == "not_sig") > 0) {
        plot_dend <- plot_dend +
            geom_segment(data=notsig_segs,
                         aes(x=x, y=y, xend=xend, yend=yend, color="not_sig"), 
                         size=1)
    }
    
    plot_dend <- plot_dend +
        scale_y_continuous(name="linkage", expand=c(.25, 0),
                           breaks=seq(0, ax_x_top, by=10^ax_x_scale)) +
                               scale_x_continuous(name="", expand=c(.1, 0),
                                                  breaks=c(),
                                                  labels=c())
    
    
    ##attach appropriate labels for colored branches along dendrogram
    plot_dend + scale_color_manual(name='Nodes',
                                   values=c('sig'='#FF1E66',
                                       'not_sig'='orange',
                                       'n_small'='#53A60A',
                                       'no_test'='#096272'),
                                   breaks=c('sig',
                                       'not_sig',
                                       'n_small',
                                       'no_test'),
                                   labels=c('significant',
                                       'not significant',
                                       'cluster too small',
                                       'untested by FWER'),
                                   drop=TRUE)
}



## #############################################################################
## #############################################################################
## helper functions

##pull label information for when hang > 0
## necessary since ggdendro package won't provide this output
.lab_height <- function(tree, heights = c()) {
    for (k in seq(length(tree))) {
        if (is.leaf(tree[[k]])) {
            heights <- c(heights, attr(tree[[k]], "height"))
        } else {
            heights <- .lab_height(tree[[k]], heights)
        }
    }
    heights
}






