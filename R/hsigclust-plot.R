#' @title plot hsigclust object
#' 
#' @description 
#' visualize results of HSigClust analysis as an annotated dendrogram
#' object using the ggdendro and ggplot2 packages
#' 
#' @details
#' some details
#' 
#' @import ggplot2 ggdendro dplyr
#' @name hsigclust-plot
#' @export
#' @author Patrick Kimes


setMethod("plot", signature(x="hsigclust", y="missing"),
          function(x, y, arg="all", ...) {
            .plot.hsigclust(x, arg, ...)
          })

.plot.hsigclust <- function(hsigclust, arg="all", ...) {
  
  hcd <- as.dendrogram(hsigclust@hc)
  
  labelColors <- c("#CDB380", "#036564", "#EB6841", "#EDC951")
  sigColor <- "#FF1E66"
  nullColor <- "gray"
  
  #using ggdendro package
  hcdata <- ggdendro::dendro_data.hclust(hsigclust@hc)
  hc_segs <- ggdendro::segment(hcdata)
  hc_labs <- ggdendro::label(hcdata)
  sig_spots <- (hsigclust@mpvalnorm[, 1] < 0.4) #need to correct, use Meinshuasen correction
  sig_linkvals <- as.factor(hsigclust@hc$height[sig_spots])
  sig_segs <- filter(hc_segs, as.factor(y) %in% sig_linkvals)
  sig_segtops <- filter(sig_segs, y == yend & x < xend)
  sig_segtops <- cbind(sig_segtops, 
                       "pval"= format(hsigclust@mpvalnorm[sig_spots, 1], 
                                      digits=3, scientific=TRUE))
  
  axis_xref <- max(hsigclust@hc$height)
  axis_xtop <- max(axis_xref)*1.25
  axis_xbot <- -max(axis_xref)*1.25/2
  axis_xscale <- floor(log10(axis_xtop))
  
  plot_dend <- ggplot() + 
    geom_segment(data=hc_segs, aes(x=x, y=y, xend=xend, yend=yend), color=nullColor) +
    geom_text(data=hc_labs, 
              aes(x=x, y=y-10, label=label, hjust=1, vjust=.5, angle=90), 
              size=3) +
    scale_y_continuous(name="linkage",
                       limits=c(axis_xbot, axis_xtop), 
                       breaks=seq(0, axis_xtop, by=10^axis_xscale)) +
    scale_x_continuous(name="",
                       breaks=c(),
                       labels=c()) +
    ggtitle("showing all p-values < 0.4") +
    theme_bw()
  
  if (length(sig_spots) > 0) {
    plot_dend <- plot_dend +
    geom_segment(data=sig_segs, aes(x=x, y=y, xend=xend, yend=yend), color=sigColor, size=1) +
      geom_text(data=sig_segtops, aes(x=x, y=y, label=pval, hjust=-0.2, vjust=-0.5),
                col=sigColor, size=4)      
  }

  plot_dend
}

