sigclust2 [![Build Status](https://travis-ci.org/pkimes/sigclust2.svg?branch=master)](https://travis-ci.org/pkimes/sigclust2)
=======================

### Contents
1. [Introduction](#intro)
2. [Tutorial](#tutorial)
    1. [Test](#test)
    2. [Plot](#plot)
3. [Details](#details)
    1. [Testing Parameters](#testparams)
    1. [Plotting Parameters](#plotparams)
4. [Miscellanea](#misc)
    1. [Installing `Rclusterpp` on OSX](#rclusterpp)
    2. [Caveats](#)
5. [References](#refs)


### <a name="intro"></a> Introduction
This package may be used to assess statistical significance in 
hierarchical clustering. To assess significance in high-dimensional data,
the approach assumes that a cluster may be well approximated by a single
Gaussian (normal) distribution. Then, given the results of hierarchical clustering,
the approach sequentially tests from the root node whether the data at each split/join
correspond to one or more Gaussian distributions. The hypothesis test performed at 
each node is based on a Monte Carlo simulation procedure, and the family-wise error
rate (FWER) is controlled across the dendrogram using a sequential testing procedure.  

A more detailed description of the approach may be found in:
Kimes PK, Liu Y, Hayes DN, Marron JS. "Statistical significance for hierarchical clustering."
(arXiv pre-print soon.)



### <a name="tutorial"></a> Tutorial
We illustrate the basic functionality of the package using a toy example. More detail on 


#### <a name="test"></a> Test


#### <a name="plot"></a> Plot



### <a name="detauls"></a> Details

#### <a name="testparams"></a> Testing Parameters

#### <a name="plotparams"></a> Plotting Parameters



### <a name="misc"></a> Miscellanea

#### <a name="rclusterpp"></a> Installing `Rclusterpp` on OSX

As described in the [`Rclusterpp` wiki][rcpp], to make use of the package's multi-threading
capabilities, a separate compiler (e.g. `gcc-4.9`)  must be installed and used to build the package. 
This essentially amounts to:

1. downloading a local `gcc` compiler (e.g. using `homebrew`)  
2. modifying your `~/.R/Makevars` file to include the following lines:
    ```{sh}
    CFLAGS += -std=c11
    CXXFLAGS += -std=c++11

    VER=-4.9
    CC=gcc$(VER)
    CXX=g++$(VER)
    SHLIB_CXXLD=g++$(VER)
    ```

3. rebuilding `Rclusterpp` and associated dependencies:
    ```{Rconsole}
    R> install.packages("Matrix")
    R> install.packages(c("Rcpp", "RcppEigen", "Rclusterpp"), type="source")
    ```


### <a name="refs"></a> References

* Kimes PK, Liu Y, Hayes DN, and Marron JS. "Statistical significance 
for hierarchical clustering." _to be submitted._
* Huang H, Liu Y, Yuan M, and Marron JS. (2014). "Statistical significance of 
clustering using soft thresholding." _Journal of Computational and Graphical Statistics_, pre-print.
* Liu Y, Hayes DN, Nobel A, and Marron JS. (2008). "Statistical significance of 
clustering for high-dimension, low–sample size data." 
_Journal of the American Statistical Association_, 103(483).


[rcpp]: https://github.com/nolanlab/Rclusterpp/wiki/Getting-Started
