\docType{methods}
\name{hsigclust-getters}
\alias{asimcindex}
\alias{asimcindex,}
\alias{clusterList}
\alias{clusterList,}
\alias{controlFWER,}
\alias{hc}
\alias{hc,}
\alias{hsigclust-getters}
\alias{icovest}
\alias{icovest,}
\alias{meigval}
\alias{meigval,}
\alias{mpval}
\alias{mpval,}
\alias{mpvalnorm}
\alias{msimeigval}
\alias{msimeigval,}
\alias{nsim}
\alias{nsim,}
\alias{plot,}
\alias{print,}
\alias{raw.data,}
\alias{show,}
\alias{summary,}
\alias{vsimbackvar}
\alias{vsimbackvar,}
\alias{xmcindex}
\title{getters for the hsigclust class}
\usage{
raw.data(x)

meigval(x)

msimeigval(x)

vsimbackvar(x)

icovest(x)

nsim(x)

asimcindex(x)

mpval(x)

mpvalnorm(x)

xmcindex(x)

clusterList(x)

hc(x)
}
\value{
raw.data: the original data used in the analysis

meigval: matrix of sample eigenvalues from observed
subtrees

msimeigval: matrix of eigenvalues used for simulating null
subtrees

vsimbackvar: vector of background noise estimates for
simulating null subtrees

icovest: the covariance estimation procedure used for the
analaysis

nsim: number of simulated null datasets for each test

asimcindex:

mpval: matrix of empirical p-values with ncols = number of
cluster indicies tested (currently 1)

mpvalnorm: matrix of empirical p-values with ncols = number
of cluster indicies tested (currently 1)

xmcindex:

clusterList:

hc:
}
\description{
Defining getters (no setters) for the hsigclust class.
}
\author{
Patrick Kimes
}

