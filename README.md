sigclust2 
=======================

### Contents
1. [Status](#status)
3. [Introduction](#intro)
4. [Example](#example)
5. [References](#refs)


### <a name="status"></a> Status
The package implements Monte Carlo based significance testing procedures for
clustering. Currently, the package includes a test for hierarchical clustering,
statisical Significance for Hierarchical Clustering  (`shc`), with the intention of
implementing a multi-cluster testing procedure. In the future, the package may be
merged with the `sigclust` package currently available through CRAN.

A short to-do list for the near future:
* `shc-class` methods:
  * `diagnostic`: implement
  * `plot`: improve flexibility (e.g. add option to print p-vals, title of plot, etc.)
* implement multicluster approach
* compile as package


### <a name="intro"></a> Introduction
This package may be used to assess the statistical significance in
hierarchical clustering. Given the results of hierarchical clustering,
the approach sequentially tests starting from the root node whether each 
split/join corresponds to "true" clustering. The hypothesis test performed at 
each node is based on the approach described in Liu et al. (2008) with 
appropriate modifications for hierarchical clustering.


<!--
### <a name="example"></a> Example
Consider the `mtcars` dataset. The HSigClust testing procedure may be 
implemented for a specific clustering procedure, e.g. euclidean dissimilarity 
and average linkage, using the call:


```r
library(sigclust2)
```

```
## Error: there is no package called 'sigclust2'
```

```r
##run HSigClust on toy dataset using Ward linkage
our_hsc <- shc(mtcars, metric="euclidean", linkage="ward")
```

```
## Error: could not find function "shc"
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
## Error: could not find function "shc"
```

We can access the p-values at each node by calling the getter function, 
`p_norm`:

```r
##only print p-values for the last 5 merges
tail(p_norm(short_hsc))
```

```
## Error: could not find function "p_norm"
```

The order of the p-values is according to the height of each branch, i.e. 
`p_norm(our_hsc[31, ])` corresponds to the highest, (n-1)st branch, at the 
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
plot(our_hsc, colGroups=mylabs, textLabs=TRUE, FWER=FALSE, alpha=1)
```

```
## Error: object 'our_hsc' not found
```

Note that the call to `plot()` returns a `ggplot` object. Therefore, we can 
easily adjust the plot using any function from the `ggplot2` package. 

-->

### <a name="refs"></a> References

* Liu Y, Hayes DN, Nobel A, and Marron JS. (2008). "Statistical significance of 
clustering for high-dimension, low–sample size data." 
_Journal of the American Statistical Association_, 103(483).
* Huang H, Liu Y, Yuan M, and Marron JS. (2014). "Statistical significance of 
clustering using soft thresholding." _Journal of Computational and Graphical Statistics_, pre-print.
* Kimes PK, Hayes DN, Liu Y, and Marron JS. "Statistical significance 
for hierarchical clustering." _submitting soon!_

