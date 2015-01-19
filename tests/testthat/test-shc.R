context("shc constructor method")


test_that("shc constructor generates shc object for appropriate input", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    ddf <- as.data.frame(dm)
    df_error <- "x must be a matrix; use as.matrix if necessary"
    
    expect_is(shc(dm), "shc")
    expect_error(shc(ddf), df_error)
})


test_that("shc metric/linkage/l parameters work properly", {
    data <- matrix(rnorm(200), ncol=2, nrow=100)
    cluster_vec <- list(
        list(metric = "euclidean", linkage = "ward.D2", l = 2, out = "pass"),
        list(metric = "euclidean", linkage = "ward", l = 2, out = "pass"),
        list(metric = "euclidean", linkage = "average", l = 2, out = "pass"),
        list(metric = "euclidean", linkage = "complete", l = 2, out = "pass"),
        list(metric = "euclidean", linkage = "single", l = 2, out = "pass"),

        list(metric = "manhattan", linkage = "average", l = 1, out = "pass"),
        list(metric = "manhattan", linkage = "single", l = 2, out = "pass"),
        list(metric = "manhattan", linkage = "complete", l = 2, out = "pass"),

        list(metric = "cor", linkage = "average", l = 2, out = "pass"),
        list(metric = "cor", linkage = "complete", l = 2, out = "pass"),
        list(metric = "cor", linkage = "single", l = 2, out = "pass")
    )
    
    for (params in cluster_vec) {
        if (params$out == "pass") {
            out <- shc(dm, metric = params$metric, linkage = params$linkage,
                       l = params$l, n_sim=5, n_min=30)

            expect_is(out, "shc")
            expect_match(out$hc_dat$method, params$linkage)
            if (params$metric != "cor") {
                expect_match(out$hc_dat$dist.method, params$metric)
            }
        }
    }
})


test_that("shc constructer alpha parameter works properly", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    alpha_error <- "invalid choice for alpha; alpha must be 0 < alpha < 1"    

    out <- shc(dm, alpha = 0)
    expect_equal(sum(out$nd_type == "sig"), 0)

    out <- shc(dm, alpha = 0.1)
    expect_is(out, "shc")

    out <- shc(dm, alpha = 1)
    expect_is(out, "shc")
    expect_equal(sum(out$nd_type == "not_sig"), 0)

    expect_error(shc(dm, alpha = 1.2), alpha_error)
})


test_that("shc constructer n_sim parameter works properly", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    n_sim_vec <- c(5, 10, 100)

    for (n_sim in n_sim_vec) {
        out <- shc(dm, n_sim = n_sim)
        expect_is(out, "shc")
        expect_equal(dim(out$ci_sim)[2], n_sim)
    }
})


test_that("shc constructer n_min parameter works properly", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    expect_is(shc(dm, n_min = 100), "shc")
    expect_is(shc(dm, n_min = 3), "shc")
    expect_error(shc(dm, n_min = 101), "n_min must be <= n")
    expect_error(shc(dm, n_min = 1), "n_min must be >= 3")
    expect_error(shc(dm, n_min = 2), "n_min must be >= 3")
})


test_that("shc constructer icovest parameter works properly", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    icovest_vec <- 1:4
    icov_warning <- "icovest should be 1, 2 or 3. Using default value: 1."
    
    expect_is(shc(dm, icovest = 1), "shc")
    expect_is(shc(dm, icovest = 2), "shc")
    expect_is(shc(dm, icovest = 3), "shc")
    expect_warning(shc(dm, icovest = 4), icov_warning)
})


test_that("shc constructer ci/null_alg parameters work properly", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    mismatch_error <- "ci = 'linkage', null_alg = '2means' cannot be specified"
    length_error <- "ci and null_alg must be of same length"
    
    expect_is(shc(dm, ci = "2CI", null_alg = "hclust"), "shc")
    expect_is(shc(dm, ci = "2CI", null_alg = "2means"), "shc")
    expect_is(shc(dm, ci = "linkage", null_alg = "hclust"), "shc")
    expect_error(shc(dm, ci = "linkage", null_alg = "2means"), mismatch_error)
    expect_error(shc(dm, ci = c("linkage", "2CI"), null_alg = "2means"), length_error)

    mout <- shc(dm, ci = c("linkage", "2CI"),
                null_alg = c("hclust", "hclust"), n_sim = 10)
    expect_equal(dim(mout$ci_sim), c(99, 10, 2))
    expect_equal(dim(mout$ci_dat), c(99, 2))
})


test_that("shc constructer ci_idx/ci_emp parameters work properly", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    ci_vec <- c("2CI", "2CI")
    null_vec <- c("hclust", "2means")
    long_error <- "invalid choice for ci_idx; ci_idx must be < length(ci)"

    expect_is(shc(dm, ci = ci_vec, null_alg = null_vec, ci_idx = 1, ci_emp = FALSE), "shc")
    expect_is(shc(dm, ci = ci_vec, null_alg = null_vec, ci_idx = 1, ci_emp = TRUE), "shc")
    expect_is(shc(dm, ci = ci_vec, null_alg = null_vec, ci_idx = 1, ci_emp = FALSE), "shc")
    expect_is(shc(dm, ci = ci_vec, null_alg = null_vec, ci_idx = 1, ci_emp = TRUE), "shc")

    expect_is(shc(dm, ci = ci_vec, null_alg = null_vec, ci_idx = 1, ci_emp = TRUE), "shc")
})


## expect_equal(str_length("a"), 1)
## expect_output(str(a), "List of 2")
## expect_match(string, "Testing")
## expect_message(library(mgcv), "This is mgcv")
## expect_warning(log(-1))
## expect_error(1 / "a")
## expect_is(lm(mpg ~ wt, data = mtcars), "lm")

