#' @title Plot hsigclust object
#'
#' @name hsigclust-plot
#' 
#' @description Visualize the results of HSigClust analysis as an annotated 
#'              dendrogram with significant branches highlighted.
#' 
#' @details This function makes use of dendrogram plotting functions made
#'          available through the ggdendro package which provides a ggplot2-like
#'          grammer for working with dendrograms.
#' @param x a hsigclust object to plot produced by a call to \code{HSCtest()}
#' @param colGroups a vector specifying group labels for the clustered objects.
#'        The vector should be in the same order as the rows of the original 
#'        data matrix. If specified, color blocks will be placed along the 
#'        bottom of the dendrogram. Useful when the samples have a priori known
#'        grouping behavior, default is \code{NULL}.
#' @param textLabs a boolean specifyin whether rowlabels should be added as 
#'        text along the bottom of the dendrogram, default is \code{TRUE}.
#' @param FWER a boolean specifying whether the FWER control procedure of 
#'        Meinshausen et al. 2010 should be used, default is \code{TRUE}. 
#'        NOTE: only has effect if \code{alphaStop} was not specified, or was
#'        set to the default value of 1 when calling \code{HSCtest()}.
#' @param alpha a double between 0 and 1 specifying the significance cutoff. If 
#'        FWER is TRUE, the FWER of the entire dendrogram is controlled at 
#'        alpha, else, each branch is tested at \code{alpha}, default is 0.05.
#'        NOTE: only has effect if \code{alphaStop} was not specified, or was
#'        set to the default value of 1, when calling \code{HSCtest()}.
#' @param hang a double value corresponding to the \code{hang} parameter for 
#'        the typical call to \code{plot} for an object of class
#'        \code{hsigclust}, default is -1.
#' 
#' @import ggplot2 ggdendro dplyr
#' @export 
#' @author Patrick Kimes


setMethod("plot", signature(x="hsigclust", y="missing"),
          function(x, y, ...) {#colLabs=TRUE, textLabs=FALSE, FWER=TRUE, alpha=0.05,
            .plot.hsigclust(x, ...)# colLabs, textLabs, FWER, alpha, ...) 
          })

.plot.hsigclust <- function(hsigclust, colGroups=NULL, textLabs=TRUE, 
                            FWER=TRUE, alpha=0.05, hang=-1) {

  
  #check to see if alphaStop was specified in original call
  if (hsigclust@alphaStop < 1) {
    cat('HSigClust procedure was applied with FWER control.')
    cat('Using original cutoff values.')
    FWER <- TRUE
    alpha <- hsigclust@alphaStop
  }
  
  #make sure significance cutoff is reasonable
  if (alpha < 0 || alpha > 1) { 
    cat('alpha must be between (0,1). Using default value of 0.05.')
    alpha <- 0.05
  }

  #size of each cluster
  #p-value cutoffs for Meinshausen correction
  n <- nrow(hsigclust@mpvalnorm)+1
  if (FWER) {
    clusterSizes <- apply(hsigclust@clusterList, 1, 
                          function(x) { length(unlist(x)) })
    cutoff <- alpha * clusterSizes/n    
  } else {
    cutoff <- rep(alpha, n-1)
  }
  
  
  allPDpairs <- rbind(cbind(hsigclust@hc$merge[,1], 1:(n-1)), 
                      cbind(hsigclust@hc$merge[,2], 1:(n-1)))
  PDmap <- data.frame(allPDpairs[allPDpairs[, 1]>0, ])
  names(PDmap) <- c("daughter", "parent")
  PDmap <- PDmap[order(PDmap$daughter), 2] #the parent of each daughter
  PDmap <- c(PDmap, n) #add final node without a parent
  
  #compute Meinshausen cutoffs for significance at alpha
  #   clusterSizes <- apply(xclust$clusterList, 1, 
  #                         function(x) { length(unlist(x)) })
  #   cutoff <- alphaStop * clusterSizes/n    
  #   
  #   treesig <- rep(TRUE, n)
  
  #### 1. NOT QUITE: STILL NEED TO CHECK PARENT--CHILD INHERITANCE
  #### 2. ALSO NEED TO ADD OPTION TO PLOT CUTOFF VALUE ABOVE/BELOW P-VALUE
  
  #colors to be used in figure
  sigColor <- "#FF1E66"
  sigtextColor <- "#A60A3C"
  nullColor <- "gray"
  nsmallColor <- "#53A60A"
  
  
  
  #using ggdendro package
  hcd <- as.dendrogram(hsigclust@hc, hang=hang)
  hcdata <- ggdendro::dendro_data(hcd)
  hc_segs <- ggdendro::segment(hcdata) #getters
  hc_labs <- ggdendro::label(hcdata) #getters
  #change label hight to be correct
  hc_labs$y <- .getLabHeight(hcd)
  
  if (!is.null(colGroups) && length(colGroups)==n) {
    hc_labs$clusters <- colGroups[hsigclust@hc$order]
  }
  
  #determine significant branches/nodes
  # -- have to determine which branches to color significant by 
  #    their linkage values, not sure if there is a better way
  sig_spots <- (hsigclust@mpvalnorm[, 1] < cutoff)
  sig_linkvals <- as.factor(hsigclust@hc$height[sig_spots])
  sig_segs <- filter(hc_segs, as.factor(y) %in% sig_linkvals)
  sig_segtops <- filter(sig_segs, (y == yend) & (x < xend))
  sig_segtops <- sig_segtops[order(sig_segtops[, 2]), ]
  sig_segtops <- cbind(sig_segtops, 
                       "pval"= format(hsigclust@mpvalnorm[sig_spots, 1], 
                                      digits=3, scientific=TRUE))
  
  #determine branches/nodes w/ not enough samples (pval=2)
  skip_spots <- (hsigclust@mpvalnorm[, 1] == 2)
  skip_linkvals <- as.factor(hsigclust@hc$height[skip_spots])
  skip_segs <- filter(hc_segs, as.factor(y) %in% skip_linkvals)
  
  
  #determine branches/nodes not tested (pval=47)

  
  axis_xref <- max(hsigclust@hc$height)
  axis_xtop <- max(axis_xref)*1.25
  axis_xbot <- -max(axis_xref)/4
  axis_xscale <- floor(log10(axis_xtop))
  y_range <- max(hc_segs$y) - min(hc_segs$y)  
  
  plot_dend <- ggplot() + 
    geom_segment(data=hc_segs, 
                 aes(x=x, y=y, xend=xend, yend=yend), 
                 color=nullColor) +
    geom_segment(data=skip_segs,
                 aes(x=x, y=y, xend=xend, yend=yend), 
                 color=nsmallColor, size=1) +
    theme_bw()

  #add appropriate title
  plot_dend <- plot_dend + 
                ggtitle(paste("showing all p-values below", 
                              alpha, "cutoff", 
                              ifelse(FWER, "(FWER corrected)","")))
  
  
  #add color group labels for objects if specified
  if (!is.null(colGroups)) {
    hc_labs$colory <- hc_labs$y - y_range/30
    plot_dend <- plot_dend +
      geom_tile(data=hc_labs,
              aes(x=x, y=colory, fill=as.factor(clusters), vjust=0), 
              alpha=1, height=y_range/30) +
      scale_fill_discrete('Labels')
  }
  
  #add text labels for each object if desired
  if (textLabs) {
    hc_labs$texty <- hc_labs$y - (2-is.null(colGroups))*y_range/30
    plot_dend <- plot_dend +
          geom_text(data=hc_labs, 
                    aes(x=x, y=texty, label=label, 
                        hjust=1, vjust=.5, angle=90), 
                    size=3)
  }
  
  
  #add colored segments to plot if any branches were significant
  if (sum(sig_spots) > 0) {
    plot_dend <- plot_dend +
      geom_segment(data=sig_segs,
                   aes(x=x, y=y, xend=xend, yend=yend), 
                   color=sigColor, size=1) +
      geom_text(data=sig_segtops, 
                aes(x=x, y=y, label=pval, hjust=-0.2, vjust=-0.5),
                col=sigtextColor, size=4)      
  }
  
  plot_dend <- plot_dend +
    scale_y_continuous(name="linkage", expand=c(.25, 0),
                       breaks=seq(0, axis_xtop, by=10^axis_xscale)) +
    scale_x_continuous(name="", expand=c(.1, 0),
                       breaks=c(),
                       labels=c())
    

  plot_dend
}



################################################################################
################################################################################
#function to pull label height/information for when hang>0
# - necessary since ggdendro package won't provide this output
.getLabHeight <- function(tree, heights=c()) {
  for (k in seq(length(tree))) {
    if (is.leaf(tree[[k]])) {
      heights <- c(heights, attr(tree[[k]], "height"))
    } else {
      heights <- getLabHeight(tree[[k]], heights)
    }
  }
  heights
}
