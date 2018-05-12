context("shc constructor method")


test_that("shc accepts only matrix data input", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## check constructor works on original data
    expect_is(shc(dm), "shc")

    ## check constructor fails with data.frame 
    err_df <- "x must be a matrix; use as.matrix if necessary"
    expect_error(shc(data.frame(dm)), err_df)
})


test_that("shc accepts common clustering algorithm parameters", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## list of valid parameter settings to test
    valid_params <- list(
        list(metric="euclidean", linkage="ward.D2"),
        list(metric="euclidean", linkage="ward"),
        list(metric="euclidean", linkage="average"),
        list(metric="euclidean", linkage="complete"),
        list(metric="euclidean", linkage="single"),
        list(metric="manhattan", linkage="average"),
        list(metric="cor", linkage="average")
    )

    ## check common hclust methods are accepted
    for (p in valid_params) {
        out <- shc(dm, n_sim=5, n_min=30, metric=p$metric, linkage=p$linkage)
        if (p$metric == "cor") {
            hc <- hclust(as.dist(1 - cor(t(dm))), method=p$linkage)
        } else {
            hc <- hclust(dist(dm, method=p$metric), method=p$linkage)
        }
        
        ## check class and clustering results
        expect_is(out, "shc")
        expect_equal(out$hc_dat$height, hc$height)
    }

    ## check minkowski distance separately
    out_d <- shc(dm, n_sim=5, n_min=30, metric="minkowski", linkage="average")
    out_3 <- shc(dm, n_sim=5, n_min=30, metric="minkowski", linkage="average", l=3)
    hc_d <- hclust(dist(dm, method="minkowski"), method="average")
    hc_3 <- hclust(dist(dm, method="minkowski", p=3), method="average")

    expect_is(out_d, "shc")
    expect_equal(out_d$hc_dat$height, hc_d$height)
    expect_is(out_3, "shc")
    expect_equal(out_3$hc_dat$height, hc_3$height)
})


test_that("shc accepts alpha values between 0 and 1", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ##  check still valid, but nothing significant
    out <- shc(dm, alpha = 0)
    expect_equal(sum(out$nd_type == "sig"), 0)

    ##  check still valid, but only 1 test
    out <- shc(dm, alpha = 1e-200)
    expect_equal(sum(out$nd_type == "not_sig"), 1)

    ## check still valid, but nothing tested
    out <- shc(dm, alpha = 1)
    expect_equal(sum(out$nd_type == "not_sig"), 0)

    ## check returns errors
    err_alpha <- "invalid choice for alpha; alpha must be 0 < alpha < 1"    
    expect_error(shc(dm, alpha = 1.2), err_alpha)
    expect_error(shc(dm, alpha = -.2), err_alpha)
})


test_that("shc accepts n_sim values", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## some simulation sizes to test
    n_sim_vec <- c(1, 5, 10, 100)

    ## check that n_sim value matches num CIs simulated
    for (n_sim in n_sim_vec) {
        out <- shc(dm, n_sim = n_sim)
        expect_equal(dim(out$ci_sim)[2], n_sim)
    }
})


test_that("shc accepts n_min values between 3 and n", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## check valid input values
    expect_is(shc(dm, n_min = 100), "shc")
    expect_is(shc(dm, n_min = 3), "shc")

    ## check invalid input values
    expect_error(shc(dm, n_min = 101), "n_min must be <= n")
    expect_error(shc(dm, n_min = 1), "n_min must be >= 3")
    expect_error(shc(dm, n_min = 2), "n_min must be >= 3")
})


test_that("shc accepts icovest values in {1, 2, 3}", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## check valid input values
    expect_is(shc(dm, n_sim=10, n_min=30, icovest = 1), "shc")
    expect_is(shc(dm, n_sim=10, n_min=30, icovest = 2), "shc")
    expect_is(shc(dm, n_sim=10, n_min=30, icovest = 3), "shc")
    
    ## check invalid input values
    warn_icov <- "icovest should be 1, 2 or 3. Using default value: 1."
    expect_warning(shc(dm, n_sim=10, n_min=30, icovest=4), warn_icov)
    expect_warning(shc(dm, n_sim=10, n_min=30, icovest=2.5), warn_icov)
    expect_warning(shc(dm, n_sim=10, n_min=30, icovest=-1), warn_icov)
})


test_that("shc accepts expected ci/null_alg parameter values", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## check that common pairs are accepted
    expect_is(shc(dm, ci = "2CI", null_alg = "hclust"), "shc")
    expect_is(shc(dm, ci = "2CI", null_alg = "2means"), "shc")
    expect_is(shc(dm, ci = "linkage", null_alg = "hclust"), "shc")
    
    ## check that linkage + 2means isn't accepted
    err_mismatch <- "ci = 'linkage', null_alg = '2means' cannot be specified"
    expect_error(shc(dm, ci = "linkage", null_alg = "2means"), err_mismatch)

    ## check that multiple ci/null_alg values can be passed
    mout <- shc(dm, ci = c("linkage", "2CI"), null_alg = c("hclust", "hclust"), n_sim = 10)
    expect_equal(dim(mout$ci_sim), c(99, 10, 2))
    expect_equal(dim(mout$ci_dat), c(99, 2))

    ## check that error occurs when ci/null_alg aren't same length
    err_length <- "ci and null_alg must be of same length"
    expect_error(shc(dm, ci = c("linkage", "2CI"), null_alg = "2means"), err_length)
})


test_that("shc accepts only valid ci_idx/ci_emp parameters", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ci_vec <- c("2CI", "2CI")
    null_vec <- c("hclust", "2means")

    ## check default parameters don't cause error
    obj <- shc(dm, n_sim=10, n_min=30, ci=ci_vec, null_alg=null_vec, ci_idx=1, ci_emp=FALSE)
    expect_is(obj, "shc")

    ## check non-default parameters don't cause error
    obj <- shc(dm, n_sim=10, n_min=30, ci=ci_vec, null_alg=null_vec, ci_idx=2, ci_emp=TRUE)
    expect_is(obj, "shc")

    ## check invalid choice causes errors
    err_long <- "invalid choice for ci_idx"
    expect_error(shc(dm, ci = ci_vec, null_alg = null_vec, ci_idx = 4, ci_emp = TRUE), err_long)
})


test_that("shc accepts user specified metric functions w/ vecmet and matmet", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(100), ncol=2, nrow=50)

    ## valid vecmet function and matrix extension
    vfun <- function(x, y) {1 - cor(x, y)}
    vfun_m <- function(x) { as.dist(outer(split(x, row(x)), split(x, row(x)), Vectorize(vfun))) }

    ## valid matmet function
    mfun <- function(x) { as.dist(1 - cor(t(x))) }

    ## check class and clustering results with vecmet
    shc_vfun <- shc(dm, n_sim=10, n_min=20, vecmet=vfun, linkage="average")
    hc_vfun <- hclust(vfun_m(dm), method="average")
    expect_is(shc_vfun, "shc")
    expect_equal(shc_vfun$hc_dat$height, hc_vfun$height)

    ## check class and clustering results with matmet
    shc_mfun <- shc(dm, n_sim=10, n_min=20, matmet=mfun, linkage="average")
    hc_mfun <- hclust(mfun(dm), method="average")
    expect_is(shc_mfun, "shc")
    expect_equal(shc_mfun$hc_dat$height, hc_mfun$height)
    
    ## check error if both vecmet and matmet are specified
    expect_error(shc(dm, vecmet=vfun, matmet=mfun, linkage="average"),
                 "only one of vecmet and matmet can be specified")

    ## check errors if vecmet, matmet are not functions
    expect_error(shc(dm, vecmet="some_text", linkage="average"),
                 paste("vecmet must be a function taking two vectors as input",
                       "and returning a real-valued dissimilarity"))
    expect_error(shc(dm, matmet="some_text", linkage="average"),
                 paste("matmet must be a function taking a data matrix as input",
                       "and returning an object of class dist"))
    
    ## invalid function specifications 
    vfun_w <- function(x, y) { warning("wrong function") }
    vfun_e <- function(x, y) { stop("wrong function") }

    ## check if warnings/errors are reurned if vecmet is a poorly specified
    expect_error(shc(dm, vecmet=vfun_w, linkage="average"),
                 "warning for vecmet specification: ")
    expect_error(shc(dm, vecmet=vfun_e, linkage="average"),
                 "error with vecmet specification: ")
})


test_that("shc accepts different bkgd noise calculations w/ bkgd_pca", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(100*100), ncol=100, nrow=100)

    ## run shc with/without PCA-based background noise estimation
    out_t <- shc(dm, n_sim=10, n_min=30, bkgd_pca=TRUE)
    out_f <- shc(dm, n_sim=10, n_min=30, bkgd_pca=FALSE)

    ## explicitly call background noise est subroutines
    ne_t <- null_eigval(dm, nrow(dm), ncol(dm), icovest=1, bkgd_pca=TRUE)
    ne_f <- null_eigval(dm, nrow(dm), ncol(dm), icovest=1, bkgd_pca=FALSE)

    ## check that PCA-based estimate is always less than standard approach
    expect_lt(ne_t$backvar, ne_f$backvar)

    ## check that shc returns same result as subroutine
    expect_equal(out_t$backvar[1], ne_t$backvar, tolerance=0.001)
    expect_equal(out_f$backvar[1], ne_f$backvar, tolerance=0.001)
})



test_that("shc accepts rcpp specification to call rclusterpp", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## run shc with/without PCA-based background noise estimation
    out_base <- shc(dm, n_sim=10, n_min=30, rcpp=FALSE)
    out_cpp <- shc(dm, n_sim=10, n_min=30, linkage="ward", rcpp=TRUE)

    ## check that clustering returns same result
    expect_equal(out_cpp$hc_dat$height, out_base$hc_dat$height)
    
    ## check that error returned with linakge = 'ward.D2', rcpp = TRUE
    expect_error(shc(dm, n_sim=10, n_min=30, rcpp=TRUE),
                 "for linkage when rcpp = TRUE.")
})

