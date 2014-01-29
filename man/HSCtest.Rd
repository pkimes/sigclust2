\name{HSCtest}
\alias{HSCtest}
\title{Run Statistical Significance of Hierarchical Clustering algorithm}
\usage{
HSCtest(x, metric, linkage, alpha = 0.05, square = FALSE, l = 2,
  nsim = 100, minObs = 10, icovest = 1)
}
\arguments{
  \item{x}{a dataset with p rows and n columns, with
  observations in columns.}

  \item{metric}{a string specifying the metric to be used
  in the hierarchical clustering procedure.  This must be a
  metric accepted by \code{dist}, e.g. "euclidean." If
  squared Euclidean distance (or the square of any other
  metric) is desired, set the \code{square} parameter to
  \code{TRUE}.}

  \item{linkage}{a string specifying the linkage to be used
  in the hierarchical clustering procedure.  This must be a
  linkage accepted by \code{hclust}, e.g. "ward."}

  \item{alpha}{a value between 0 and 1 specifying the
  desired level of the test. If no FWER control is desired,
  simply set alpha to 1. The default is 0.05.}

  \item{square}{a logical specifying whether to square the
  dissimilarity matrix produced by specified \code{metric}.
  This is necessary, for example, in order to implement
  Ward's minimum variance method with squared Euclidean
  metric. (Honestly can't think of any other situations
  when you'd want to square the diss. matrix.) Default is
  FALSE.}

  \item{l}{an integer value specifying the power of the
  Minkowski distance, if used, default is 2.}

  \item{nsim}{a numeric value specifying the number of
  simulations for SigClust testing.  The default is to run
  100 simulations at each merge.}

  \item{minObs}{an integer specifying the minimum number of
  observations needed to calculate a p-value, default is
  10.}

  \item{icovest}{a numeric value specifying the covariance
  estimation method: 1. Use a soft threshold method as
  constrained MLE (default); 2. Use sample covariance
  estimate (recommended when diagnostics fail); 3. Use
  original background noise thresholded estimate (from Liu
  et al., (2008)) ("hard thresholding") as described in the
  \code{sigclust} package documentation.}
}
\value{
The function returns a \code{hsigclust} object containing
the resulting p-values.  The print method call will output
a dendrogram with the corresponding p-values placed at each
merge.
}
\description{
\code{HSCtest} first computes the desired dendrogram, and
quantifies the statistical significance of clustering at
each branching along the tree using a hypothesis testing
procedure. The code is written so that various cluster
inidices can be easily introduced by writting new CI
function and adding to .initcluster(), .simcluster() and
incrementing "nCIs" (more elegant solution?). This function
WILL make use \code{Rclusterpp.hclust} package for all
choices of clustering except \code{dist="cor"}, for which
\code{WGCNA::cor} will be used with \code{stats::hclust}.
}
\details{
The function extends the \code{sigclust} idea to the
hierarchical setting by modifying the clustering procedure
employed to compute the null distribution of cluster
indices.
}
\author{
Patrick Kimes
}

