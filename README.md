sigclust2 [![Build Status](https://travis-ci.org/pkimes/sigclust2.svg?branch=master)](https://travis-ci.org/pkimes/sigclust2)
=======================

### Contents
1. [Status](#status)
2. [Introduction](#intro)
3. [Installing `Rclusterpp` on OSX](#rclusterpp)
4. [References](#refs)


### <a name="status"></a> Status
The package implements Monte Carlo based significance testing procedures for
clustering. Currently, the package includes a test for hierarchical clustering,
statisical Significance for Hierarchical Clustering  (`shc`), with the intention of
implementing a multi-cluster testing procedure. In the future, the package may be
merged with the `sigclust` package currently available through CRAN. The package was
originally written using [S4 objects](http://adv-r.had.co.nz/S4.html), but was switched
to the simpler [S3](http://adv-r.had.co.nz/S3.html) style.

A short to-do list for the near future:
* `shc-class` methods:
  * implement `diagnostic`
  * add flexibility to `plot`


### <a name="intro"></a> Introduction
This package may be used to assess the statistical significance in
hierarchical clustering. Given the results of hierarchical clustering,
the approach sequentially tests starting from the root node whether each 
split/join corresponds to "true" clustering. The hypothesis test performed at 
each node is based on the approach described in Liu et al. (2008) with 
appropriate modifications for hierarchical clustering.


### <a name="rclusterpp"></a> Installing `Rclusterpp` on OSX

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

3. rebuilding `Rclusterpp` and associated dependencies
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
