SigClust2 
=======================

### Contents
1. [Status](#status)
2. [Introduction](#intro)
3. [Example](#example)
4. [References](#refs)


### <a name="status"></a> Status
`R` implementation of various extensions to the CRAN `sigclust` package.
Currently, the package only includes a barely-working form of a hierarchical 
extension, HSigClust, with the intention of expanding the package to also 
include a multi-cluster extension, KSigClust.

A short to-do list for the near future:
* clean up hsigclust methods:
  * `show`: produce more useful output
  * `print`: copy `show` output
  * `HSCtest`: use `Rclusterpp.hclust` and `WGCNA::cor`
  * `diagnostics`: make work
  * `summary`: produce more useful output
  * `FWERcontrol`: make work
  * `plot`: update w/ `FWERcontrol`, add labeling option
* revise `hsigclust` class to be "lighter"
* translate KSCtest from Matlab to `R`
 * complete vignette - switch to `knitr`?
 * update README along the way (this!)


### <a name="intro"></a> Introduction
This package may be used to assess the statistical significance of clustering in
hieararchical clustering. Given the results of hierarchical clustering,
the approach sequentially tests starting from the root node whether each 
split/join corresponds to "true" clustering. The hypothesis test performed at each
node is based on the approach described in Liu et al. (2008) with appropriate
modifications for hierarchical clustering. The work is ongoing, and may involve
substantial changes to the current approach (and code).


### <a name="example"></a> Example
Consider the `mtcars` dataset. The HSigClust testing procedure may be implemented  
for a specific clustering procedure, e.g. euclidean dissimilarity and average linkage, 
using the call:


```r
library(SigClust2)
our_hsc <- HSCtest(mtcars, metric = "euclidean", linkage = "average")
```



A quick way to check the results is to simply `plot()` the output. The corresponding 
dendrogram is returned with significant splits appropriately labeled:


```r
plot(our_hsc)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 



### <a name="refs"></a> References

* Liu Y, Hayes DN, Nobel A, and Marron JS. (2008). "Statistical significance of clustering for high-dimension, low–sample size data." _Journal of the American Statistical Association_, 103(483).
* Huang H, Liu Y, Yuan M, and Marron JS. (2013). "Statistical significance of clustering using soft thresholding." _arXiv preprint [arXiv:1305.5879]_.
* Kimes P, Hayes DN, Liu Y, and Marron JS. "HSigClust: Statistical significance of hierarchical clustering." _In preparation_.

[arXiv:1305.5879]: http://arxiv.org/abs/1305.5879
