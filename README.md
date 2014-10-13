sigclust2 [![Build Status](https://travis-ci.org/pkimes/sigclust2.svg)](https://travis-ci.org/pkimes/sigclust2)
=======================

### Contents
1. [Status](#status)
3. [Introduction](#intro)
4. [Example](#example)
5. [Installing `Rclusterpp` on OSX](#rclusterpp)
6. [References](#refs)


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


<!--
### <a name="example"></a> Example
Consider the `mtcars` dataset. The SHC testing procedure may be 
implemented for a specific clustering procedure, e.g. euclidean dissimilarity 
and average linkage, using the call:


```r
library(sigclust2)
##run HSigClust on toy dataset using Ward linkage
our_hsc <- shc(mtcars, metric="euclidean", linkage="ward")
```

```
## Error: !!  x must be matrix, use as.matrix(x) if necessary.  !!
```

In the above call to `HSCtest()` we use the default value of `alpha = 1`
which results in the procedure testing at all branches along the dendrogram.
Alternatively, we may have specified `alpha = 0.05` for the testing procedure
to iteratively test from the top using a FWER control stopping procedure 
originally described in Meinshausen et al. 2010. 


```r
##run SHC on toy dataset with FWER control at 0.05
## algorithm will skip all tests ignored by sequential
## FWER procedure.
short_hsc <- shc(mtcars, metric="euclidean", linkage="ward", 
                 alpha=0.05)
```

```
## Error: !!  x must be matrix, use as.matrix(x) if necessary.  !!
```

We can access the p-values at each node by calling the getter function, 
`p_norm`:

```r
##only print p-values for the last 5 merges
tail(short_hsc$p_norm)
```

```
## Error: object 'short_hsc' not found
```

The order of the p-values is according to the height of each branch, i.e. 
`our_hsc$p_norm[31, ]` corresponds to the highest, (n-1)st branch, at the 
top of the dendrogram. p-values of `2` correspond to branches not having enough
samples to test according to the `min_n` parameter. p-values of `-1`
correspond to branches skipped according to the FWER control procedure (these 
will supercede `2` values).  

A quick way to check the results is to simply `plot` the output. The 
corresponding dendrogram is returned with significant splits appropriately 
labeled:


```r
plot(our_hsc)
```

```
## Error: object 'our_hsc' not found
```

Other plotting options are possible. Suppose we are interested in looking at 
how Mercedes cars might be distributed along the dendrogram.


```r
##extract car maker names
makers <- sapply(strsplit(rownames(mtcars), ' '), '[[', 1)
mylabs <- ifelse(makers == "Merc",
                 "Mercedes", "other")

##plot dendrogram showing all p-values and include "mercedes" label
plot(our_hsc, groups=mylabs, labs=TRUE, fwer=FALSE, alpha=1)
```

```
## Error: object 'our_hsc' not found
```

Note that the call to `plot()` returns a `ggplot` object. Therefore, we can 
easily adjust the plot using any function from the `ggplot2` package. 

-->


### <a name="rclusterpp"></a> Installing `Rclusterpp` on OSX

As described in the [`Rclusterpp` wiki](rcpp), to make use of the package's multi-threading
capabilities, a separate compiler must be installed and used to compile the package. 
This essentially ammounts to:
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
for hierarchical clustering." _submitted._
* Huang H, Liu Y, Yuan M, and Marron JS. (2014). "Statistical significance of 
clustering using soft thresholding." _Journal of Computational and Graphical Statistics_, pre-print.
* Liu Y, Hayes DN, Nobel A, and Marron JS. (2008). "Statistical significance of 
clustering for high-dimension, low–sample size data." 
_Journal of the American Statistical Association_, 103(483).



[rcpp]: https://github.com/nolanlab/Rclusterpp/wiki
