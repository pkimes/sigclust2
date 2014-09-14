\docType{methods}
\name{hsigclust-getters}
\alias{FWERcutoffs}
\alias{asimcindex}
\alias{clusterList}
\alias{controlFWER}
\alias{hc}
\alias{hsigclust-getters}
\alias{meigval}
\alias{mpval}
\alias{mpvalnorm}
\alias{msimeigval}
\alias{raw.data}
\alias{vsimbackvar}
\alias{xmcindex}
\title{getters for the hsigclust class}
\usage{
raw.data(x)

meigval(x)

msimeigval(x)

vsimbackvar(x)

asimcindex(x)

mpval(x)

mpvalnorm(x)

xmcindex(x)

clusterList(x)

hc(x)

FWERcutoffs(hsc, alpha)
}
\value{
raw.data: the original data used in the analysis

meigval: matrix of sample eigenvalues from observed
subtrees

msimeigval: matrix of eigenvalues used for simulating null
subtrees

vsimbackvar: vector of background noise estimates for
simulating null subtrees

asimcindex:

mpval: matrix of empirical p-values with ncols = number of
cluster indicies tested (currently 1)

mpvalnorm: matrix of empirical p-values with ncols = number
of cluster indicies tested (currently 1)

xmcindex:

clusterList:

hc:

FWERcutoffs:
}
\description{
Defining getters (no setters) for the hsigclust class.
}
\author{
Patrick Kimes
}

