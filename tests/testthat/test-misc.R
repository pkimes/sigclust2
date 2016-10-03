library("sigclust2")
context("miscellaneous exported functions")

test_that("fwer_cutoff called on shc object gives same results as internal function", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    alpha <- 0.10
    out <- shc(dm)
    expect_equal(fwer_cutoff(out, alpha),
                 sigclust2:::fwer_cutoff(out$idx_hc, alpha))
})


test_that("null_eigval procedure can be called externally", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    onames <- c("eigval_dat", "backvar", "eigval_sim")
    dim_error <- "Wrong size of matrix x!"
    covest_error <- "covest must be 1, 2 or 3"
    
    expect_is(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = TRUE), "list")
    expect_is(null_eigval(dm, 100, 2, icovest = 2, bkgd_pca = TRUE), "list")
    expect_is(null_eigval(dm, 100, 2, icovest = 3, bkgd_pca = TRUE), "list")
    expect_is(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = FALSE), "list")
    expect_equal(names(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = FALSE)), onames)
    expect_error(null_eigval(dm, 2, 100, bkgd_pca = FALSE), dim_error)
    expect_error(null_eigval(dm, 100, 2, icovest = 0, bkgd_pca = FALSE), covest_error)

    expect_equal(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = FALSE)$eigval_sim,
                 null_eigval(dm, 100, 2, icovest = 3, bkgd_pca = FALSE)$eigval_sim)
})


