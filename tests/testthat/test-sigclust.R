context("sigclust constructor method")


test_that("sigclust accepts only matrix data input", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## check constructor works on original data
    expect_is(sigclust(dm), "sigclust")

    ## check constructor fails with data.frame 
    df_error <- "x must be a matrix; use as.matrix if necessary"
    expect_error(sigclust(data.frame(dm)), df_error)
})


test_that("sigclust accepts n_sim values", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## some simulation sizes to test
    n_sim_vec <- c(1, 5, 10, 100)

    ## check that n_sim value matches num CIs simulated
    for (n_sim in n_sim_vec) {
        out <- sigclust(dm, n_sim = n_sim)
        expect_equal(length(out$ci_sim), n_sim)
    }
})


test_that("sigclust accepts icovest values in {1, 2, 3}", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## check valid input values
    expect_is(sigclust(dm, n_sim=10, icovest = 1), "sigclust")
    expect_is(sigclust(dm, n_sim=10, icovest = 2), "sigclust")
    expect_is(sigclust(dm, n_sim=10, icovest = 3), "sigclust")
    
    ## check invalid input values
    warn_icov <- "icovest should be 1, 2 or 3. Using default value: 1."
    expect_warning(sigclust(dm, n_sim=10, icovest=4), warn_icov)
    expect_warning(sigclust(dm, n_sim=10, icovest=2.5), warn_icov)
    expect_warning(sigclust(dm, n_sim=10, icovest=-1), warn_icov)
})


test_that("sigclust accepts different bkgd noise calculations w/ bkgd_pca", {
    ## simple dataset
    dm <- matrix(rnorm(100*100), ncol=100, nrow=100)

    ## run sigclust with/without PCA-based background noise estimation
    out_t <- sigclust(dm, n_sim=10, bkgd_pca=TRUE)
    out_f <- sigclust(dm, n_sim=10, bkgd_pca=FALSE)

    ## explicitly call background noise est subroutines
    ne_t <- null_eigval(dm, nrow(dm), ncol(dm), icovest=1, bkgd_pca=TRUE)
    ne_f <- null_eigval(dm, nrow(dm), ncol(dm), icovest=1, bkgd_pca=FALSE)

    ## check that PCA-based estimate is always less than standard approach
    expect_lt(ne_t$backvar, ne_f$backvar)

    ## check that sigclust returns same result as subroutine
    expect_equal(out_t$backvar, ne_t$backvar)
    expect_equal(out_f$backvar, ne_f$backvar)
})


test_that("sigclust accepts user-specified cluster labels", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## valid and invalid label options
    numlabs <- rep(1:2, each=50)
    letlabs <- rep(letters[1:2], each=50)
    onelabs <- rep(1, 100)
    
    ## check that specified labels used and gen lower CI
    expect_lt(sigclust(dm, n_sim=10)$ci_dat,
              sigclust(dm, n_sim=10, labels=numlabs)$ci_dat)

    ## check that labels cannot be letters
    expect_error(sigclust(dm, n_sim=10, labels=letlabs),
                 "labels must be a n-vector of 1s and 2s")

    ## check that labels cannot just be single cluster
    expect_error(sigclust(dm, n_sim=10, labels=onelabs),
                 "labels must be a n-vector of 1s and 2s")
})



