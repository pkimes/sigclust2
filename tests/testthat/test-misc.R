context("miscellaneous exported functions")


test_that("fwer_cutoff called on shc object gives same results as internal function", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm)

    ## check fwer_cutoff class methods return same result 
    expect_equal(fwer_cutoff(out, 0.10),
                 sigclust2:::fwer_cutoff.matrix(out$idx_hc, 0.10))
})


test_that("null_eigval procedure can be called externally", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    
    ## check outputs list
    expect_is(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = TRUE), "list")
    expect_is(null_eigval(dm, 100, 2, icovest = 2, bkgd_pca = TRUE), "list")
    expect_is(null_eigval(dm, 100, 2, icovest = 3, bkgd_pca = TRUE), "list")
    expect_is(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = FALSE), "list")

    ## check output structure
    expect_equal(names(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = FALSE)),
                 c("eigval_dat", "backvar", "eigval_sim"))

    ## check errors and warnings
    error_dim <- "Wrong size of matrix x!"
    warn_covest <- "icovest should be 1, 2 or 3. Using default value: 1."
    expect_error(null_eigval(dm, 2, 100, bkgd_pca = FALSE), error_dim)
    expect_warning(null_eigval(dm, 100, 2, icovest = 0, bkgd_pca = FALSE),
                   warn_covest)
    
    ## check icovest 1,3 same for low-dim
    expect_equal(null_eigval(dm, 100, 2, icovest = 1, bkgd_pca = FALSE)$eigval_sim,
                 null_eigval(dm, 100, 2, icovest = 3, bkgd_pca = FALSE)$eigval_sim)
})


test_that("shcutree returns cluster labels for shc object", {
    set.seed(100)

    ## per cluster dimensions
    n <- 100
    p <- 2

    ## simple null dataset with 1 cluster
    dm1 <- matrix(rnorm(p*n), ncol=p, nrow=n)
    out1 <- shc(dm1, alpha=0.1,
                ci=c("2CI", "linkage"),
                null_alg=c("hclust", "hclust"))
    hcout1 <- hclust(dist(dm1, method="euclidean"), method="ward.D2")

    ## simple non-null dataset with 2 clusters
    dm2 <- rbind(matrix(rnorm(p*n, mean=-3), ncol=p, nrow=n),
                 matrix(rnorm(p*n, mean=3), ncol=p, nrow=n))
    out2 <- shc(dm2, alpha=0.1,
                ci=c("2CI", "linkage"),
                null_alg=c("hclust", "hclust"))
    hcout2 <- hclust(dist(dm2, method="euclidean"), method="ward.D2")

    ## check errors for invalid input parameters
    warn_alpha <- "invalid choice for alpha; alpha must be 0 < alpha < 1"
    warn_ciidx <- "invalid choice for ci_idx; ci_idx must be a valid column index from p_emp or p_norm"
    expect_error(shcutree(out1, alpha=1.5), warn_alpha)
    expect_error(shcutree(out1, alpha=-1), warn_alpha)
    expect_error(shcutree(out1, ci_idx=3), warn_ciidx)

    ## check warnings for params disagree w/ original shc params
    warn_newalpha <- "shc constructed using smaller alpha than specified; using alpha = 0.1"
    warn_newciidx <- "shc constructed using alpha < 1, using ci_idx from obj\\$in_args"
    warn_newciemp <- "shc constructed using alpha < 1, using ci_emp from obj\\$in_args"
    expect_warning(shcutree(out1, alpha=.5), warn_newalpha)
    expect_warning(shcutree(out1, ci_idx=2), warn_newciidx)
    expect_warning(shcutree(out1, ci_emp=TRUE), warn_newciemp)

    ## check returns vector of labels with valid parameters
    expect_is(shcutree(out1), "numeric")
    expect_length(shcutree(out1), n)
    expect_is(shcutree(out1, alpha=.05), "numeric")
    expect_length(shcutree(out1, alpha=.05), n)
    expect_is(shcutree(out1, ci_idx=1), "numeric")
    expect_length(shcutree(out1, ci_idx=1), n)

    ## check one case for 2 cluster example
    expect_is(shcutree(out2), "numeric")
    expect_length(shcutree(out2), 2*n)
    
    ## check that cuts at correct height, prints correct labels
    ri1 <- WGCNA::randIndex(table(shcutree(out1, 0.001), rep(1, n)), adjust=FALSE)
    ri2 <- WGCNA::randIndex(table(shcutree(out2, 0.001), cutree(hcout2, 2)), adjust=FALSE)
    expect_equal(ri1, 1) 
    expect_equal(ri2, 1) 
})


test_that("summary and print return expected results", {
    set.seed(100)

    ## non-default parameters for testing print/summary
    .metric <- "manhattan"
    .linkage <- "ward.D"
    .alpha <- .88
    .icovest <- 3
    .ci <- c("linkage", "linkage")
    .null_alg <- c("hclust", "2CI")
    .ci_idx <- 2
    .n_sim <- 60
    .n_min <- 15

    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    ## check that non-default params are returned with 'summary'
    out <- shc(dm,
               metric = .metric,
               linkage = .linkage,
               alpha = .alpha,
               icovest = .icovest,
               ci = .ci,
               null_alg = .null_alg,
               ci_idx = .ci_idx,
               n_sim = .n_sim,
               n_min = .n_min)

    ## check that non-default params are returned with 'print'
    expect_output(print(out), "shc object created using")
    expect_output(print(out), paste0("dissimilarity = ", .metric))
    expect_output(print(out), paste0("linkage = ", .linkage))
    expect_output(print(out), paste0("n_sim = ", .n_sim))
    expect_output(print(out), paste0("icovest = ", .icovest))
    expect_output(print(out), paste0("ci = ", paste0(.ci, collapse=", ")))
    expect_output(print(out), paste0("null_alg = ", paste0(.null_alg, collapse=", ")))
    expect_output(print(out), paste0("n_min = ", .n_min))
    expect_output(print(out), paste0("FWER control = ", (.alpha < 1)))
    expect_output(print(out), paste0("alpha = ", .alpha))
    expect_output(print(out), paste0("ci_idx = ", .ci_idx))

    ## check that non-default params are returned with 'summary'
    expect_output(summary(out), "shc object created using")
    expect_output(summary(out), paste0("dissimilarity = ", .metric))
    expect_output(summary(out), paste0("linkage = ", .linkage))
    expect_output(summary(out), paste0("n_sim = ", .n_sim))
    expect_output(summary(out), paste0("icovest = ", .icovest))
    expect_output(summary(out), paste0("ci = ", paste0(.ci, collapse=", ")))
    expect_output(summary(out), paste0("null_alg = ", paste0(.null_alg, collapse=", ")))
    expect_output(summary(out), paste0("n_min = ", .n_min))
    expect_output(summary(out), paste0("FWER control = ", (.alpha < 1)))
    expect_output(summary(out), paste0("alpha = ", .alpha))
    expect_output(summary(out), paste0("ci_idx = ", .ci_idx))
})
