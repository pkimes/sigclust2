\name{HSCtest}
\alias{HSCtest}
\title{Statistical Significance of Hierarchical Clustering (HSigClust)
       algorithm}
\usage{
hsc_cars <- HSCtest(mtcars, metric="euclidean", linkage="single")
tail(mpvalnorm(hsc_cars), 10)
}
\arguments{
  \item{x}{a dataset with n rows and p columns, with
  observations in columns.}

  \item{metric}{a string specifying the metric to be used
  in the hierarchical clustering procedure. This must be a
  metric accepted by \code{dist}, e.g. "euclidean." If
  squared Euclidean distance (or the square of any other
  metric) is desired, set the \code{square} parameter to
  \code{TRUE}.}

  \item{linkage}{a string specifying the linkage to be used
  in the hierarchical clustering procedure. This must be a
  linkage accepted by \code{hclust}, e.g. "ward."}

  \item{square}{a logical specifying whether to square the
  dissimilarity matrix produced by specified \code{metric}.
  This is necessary, for example, in order to implement
  Ward's minimum variance method with squared Euclidean
  metric. (Honestly can't think of any other situations
  when you'd want to square the diss. matrix.) Default is
  FALSE.}

  \item{l}{an integer value specifying the power of the
  Minkowski distance, if used, default is 2.}

  \item{alphaStop}{a value between 0 and 1 specifying the
  desired level of the test. If no FWER control is desired,
  simply set alpha to 1. The default is 1. The testing
  procedure will terminate when no branches meet the
  corresponding FWER control threshold. This procedure may
  substantially speed up the procedure by reducing the
  number of tests considered. NOTE: if p-values are desired
  for all branches, the FWER cutoffs may be provided a
  posteriori by a call to \code{FWERcutoffs()}.
  Additionally, they may be specified at \code{plot()}
  using the \code{alpha} and \code{FWER} parameters.}

  \item{nsim}{a numeric value specifying the number of
  simulations for SigClust testing. The default is to run
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

  \item{verb}{a logical value specifying whether the method
  should print out when testing completes along each node
  along the dendrogram, by default set to FALSE.}

  \item{testCIs}{a string vector specifying the cluster
  indices to be used for testing along the dendrogram.
  Currently, options include: "2CI", "linkage". Default is
  "2CI".}

  \item{testNulls}{a string vector specifying the
  clustering approach that should be used at each node as
  the comparison. Currently, options include: "2means",
  "hclust". Default is "hclust". Note, testNulls and
  testCIs must be of equal length.}

  \item{cutoffCI}{a single value between 1 and
  \code{length(testCIs)} specifiying which CI to use for
  the FWER stopping rule, default is 1.  This only has an
  effect if alphaStop is specified to a non-default value.}

  \item{gpu}{a logical value specifying whether a GPU
  process is available and should be used. This calls
  \code{gputools::gpuDistClust} to perform the clustering
  within each simulation loop. By default set to FALSE.}
}
\value{
The function returns a \code{hsigclust} object containing
the resulting p-values. The print method call will output a
dendrogram with the corresponding p-values placed at each
merge.
}
\description{
\code{HSCtest} first computes the desired dendrogram, and
quantifies the statistical significance of clustering at
each branching along the tree using a hypothesis testing
procedure.  The code is written so that various cluster
indices can be easily introduced by writting new CI
function and adding to .initcluster(), .simcluster() and
incrementing "nCIs" (more elegant solution?). When
possible, this function makes use of a C++ implementation
of hierarchical clustering available through the
\code{Rclusterpp.hclust} package for the case of clustering
by \code{dist="cor"}, we make use of \code{WGCNA::cor} with
the usual \code{stats::hclust}.
}
\details{
The function expands on the \code{sigclust} idea to the
hierarchical setting by modifying the clustering procedure
employed to compute the null distribution of cluster
indices.
}
\author{
Patrick Kimes
}

