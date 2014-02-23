\docType{methods}
\name{hsigclust-plot}
\alias{hsigclust-plot}
\title{Plot hsigclust object}
\usage{
\S4method{plot}{hsigclust,missing}(x, y, ...)
}
\arguments{
  \item{x}{a hsigclust object to plot produced by a call to
  \code{HSCtest()}}

  \item{colGroups}{a vector specifying group labels for the
  clustered objects.  The vector should be in the same
  order as the rows of the original data matrix. If
  specified, color blocks will be placed along the bottom
  of the dendrogram. Useful when the samples have a priori
  known grouping behavior, default is \code{NULL}.}

  \item{textLabs}{a boolean specifyin whether rowlabels
  should be added as text along the bottom of the
  dendrogram, default is \code{TRUE}.}

  \item{FWER}{a boolean specifying whether the FWER control
  procedure of Meinshausen et al. 2010 should be used,
  default is \code{TRUE}.  NOTE: only has effect if
  \code{alphaStop} was not specified, or was set to the
  default value of 1 when calling \code{HSCtest()}.}

  \item{alpha}{a double between 0 and 1 specifying the
  significance cutoff. If FWER is TRUE, the FWER of the
  entire dendrogram is controlled at alpha, else, each
  branch is tested at \code{alpha}, default is 0.05.  NOTE:
  only has effect if \code{alphaStop} was not specified, or
  was set to the default value of 1, when calling
  \code{HSCtest()}.}

  \item{hang}{a double value corresponding to the
  \code{hang} parameter for the typical call to \code{plot}
  for an object of class \code{hsigclust}, default is -1.}

  \item{y}{the y coordinates of points in the plot,
  \emph{optional} if \code{x} is an appropriate structure.}

  \item{...}{Arguments to be passed to methods, such as
  \link{graphical parameters} (see \code{\link{par}}).
  Many methods will accept the following arguments:
  \describe{ \item{\code{type}}{what type of plot should be
  drawn.  Possible types are \itemize{ \item \code{"p"} for
  \bold{p}oints, \item \code{"l"} for \bold{l}ines, \item
  \code{"b"} for \bold{b}oth, \item \code{"c"} for the
  lines part alone of \code{"b"}, \item \code{"o"} for both
  \sQuote{\bold{o}verplotted}, \item \code{"h"} for
  \sQuote{\bold{h}istogram} like (or \sQuote{high-density})
  vertical lines, \item \code{"s"} for stair \bold{s}teps,
  \item \code{"S"} for other \bold{s}teps, see
  \sQuote{Details} below, \item \code{"n"} for no plotting.
  } All other \code{type}s give a warning or an error;
  using, e.g., \code{type = "punkte"} being equivalent to
  \code{type = "p"} for S compatibility.  Note that some
  methods, e.g. \code{\link{plot.factor}}, do not accept
  this.  } \item{\code{main}}{an overall title for the
  plot: see \code{\link{title}}.} \item{\code{sub}}{a sub
  title for the plot: see \code{\link{title}}.}
  \item{\code{xlab}}{a title for the x axis: see
  \code{\link{title}}.} \item{\code{ylab}}{a title for the
  y axis: see \code{\link{title}}.} \item{\code{asp}}{the
  \eqn{y/x} aspect ratio, see \code{\link{plot.window}}.} }
  }
}
\description{
Visualize the results of HSigClust analysis as an annotated
dendrogram with significant branches highlighted.
}
\details{
This function makes use of dendrogram plotting functions
made available through the ggdendro package which provides
a ggplot2-like grammer for working with dendrograms.
}
\author{
Patrick Kimes
}

