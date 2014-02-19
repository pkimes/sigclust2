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
#' @param x a hsigclust object to plot produced by a call to HSCtest()
#' @param colGroups a vector specifying group labels for the clustered objects.
#'        The vector should be in the same order as the rows of the original 
#'        data matrix. If specified, color blocks will be placed along the 
#'        bottom of the dendrogram. Useful when the samples have a priori known
#'        grouping behavior, default is NULL.
#' @param textLabs a boolean specifyin whether rowlabels should be added as 
#'        text along the bottom of the dendrogram, default is FALSE.
#' @param FWER a boolean specifying whether the FWER control procedure of 
#'        Meinshausen et al. 2010 should be used, default is TRUE.
#' @param alpha a double between 0 and 1 specifying the significance cutoff. If 
#'        FWER is TRUE, the FWER of the entire dendrogram is controlled at 
#'        alpha, else, each branch is tested at alpha, default is 0.05.
#' 
#' @import ggplot2 ggdendro dplyr
#' @export 
#' @author Patrick Kimes


setMethod("plot", signature(x="hsigclust", y="missing"),
          function(x, y, ...) {#colLabs=TRUE, textLabs=FALSE, FWER=TRUE, alpha=0.05,
            .plot.hsigclust(x, ...)# colLabs, textLabs, FWER, alpha, ...) 
          })

.plot.hsigclust <- function(hsigclust, colGroups=NULL, textLabs=FALSE, 
                            FWER=TRUE, alpha=0.05) {

  #make sure significance cutoff is reasonable
  if (alpha < 0 || alpha > 1) { 
    cat('alpha must be between (0,1). Using default value of 0.05.')
    alpha <- 0.05 
  }

  #size of each cluster
  n <- nrow(hsigclust@mpval)+1
  
  #p-value cutoffs for Meinshausen correction
  if (FWER) {
    clusterSizes <- apply(hsigclust@clusterList, 1, 
                          function(x) { length(unlist(x)) })
    cutoff <- alpha * clusterSizes/n    
  } else {
    cutoff <- rep(alpha, n-1)
  }
  
  
  #colors to be used in figure
  sigColor <- "#FF1E66"
  nullColor <- "gray"
  
  
  #using ggdendro package
  hcdata <- ggdendro::dendro_data.hclust(hsigclust@hc)
  hc_segs <- ggdendro::segment(hcdata)
  hc_labs <- ggdendro::label(hcdata)
  if (!is.null(colGroups) && length(colGroups)==n) {
    hc_labs$clusters <- colGroups[hsigclust@hc$order]
  }
  sig_spots <- (hsigclust@mpvalnorm[, 1] < cutoff)
  sig_linkvals <- as.factor(hsigclust@hc$height[sig_spots])
  sig_segs <- filter(hc_segs, as.factor(y) %in% sig_linkvals)
  sig_segtops <- filter(sig_segs, (y == yend) & (x < xend))
  sig_segtops <- sig_segtops[order(sig_segtops[, 2]), ]
  sig_segtops <- cbind(sig_segtops, 
                       "pval"= format(hsigclust@mpvalnorm[sig_spots, 1], 
                                      digits=3, scientific=TRUE))
  
  
  axis_xref <- max(hsigclust@hc$height)
  axis_xtop <- max(axis_xref)*1.25
  axis_xbot <- -max(axis_xref)/4
  axis_xscale <- floor(log10(axis_xtop))
  y_range <- max(hc_segs$y) - min(hc_segs$y)  
  
  plot_dend <- ggplot() + 
    geom_segment(data=hc_segs, 
                 aes(x=x, y=y, xend=xend, yend=yend), 
                 color=nullColor) +
    scale_y_continuous(name="linkage",
                       limits=c(axis_xbot, axis_xtop), 
                       breaks=seq(0, axis_xtop, by=10^axis_xscale)) +
    scale_x_continuous(name="",
                       breaks=c(),
                       labels=c()) +
    ggtitle("showing all p-values below FWER cutoff") +
    theme_bw()
  
  
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
  
  if (length(sig_spots) > 0) {
    plot_dend <- plot_dend +
      geom_segment(data=sig_segs,
                   aes(x=x, y=y, xend=xend, yend=yend), 
                   color=sigColor, size=1) +
      geom_text(data=sig_segtops, 
                aes(x=x, y=y, label=pval, hjust=-0.2, vjust=-0.5),
                col=sigColor, size=4)      
  }

  plot_dend
}

