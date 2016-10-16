context("shc diagnostics plotting")


test_that("plot default call works properly", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, n_sim = 10, n_min = 30, alpha = 1)

    ## check single plots return ggplot
    expect_is(diagnostic(out, pty="background"), "ggplot")
    expect_is(diagnostic(out, pty="qq"), "ggplot")
    expect_is(diagnostic(out, pty="covest"), "ggplot")
    expect_is(diagnostic(out, pty="pvalue"), "ggplot")

    ## check 'all' plot type saves to default file name
    expect_is(diagnostic(out, pty="all"), "integer")
    expect_true(file.exists("shc_diagnostic.pdf"))
    unlink("shc_diagnostic.pdf")
    expect_false(file.exists("shc_diagnostic.pdf"))

    ## check that pty can only be single unique plot type
    err_2many <- "can only specify one plot type"
    err_ambig <- "specified plot type matches more than one plot type"
    err_none <- "pty is not a valid value, please specify one of"
    expect_error(diagnostic(out, pty=c("qq", "background")), err_2many)
    expect_error(diagnostic(out, pty="a"), err_ambig)
    expect_error(diagnostic(out, pty="bad_name"), err_none)
})


test_that("plot accepts fname input", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, n_sim = 10, n_min = 30, alpha = 1)

    ## check that pdf suffix is added to fname
    expect_is(diagnostic(out, K=1, pty="all", fname="plot_1.jpeg"), "integer")
    expect_true(file.exists("plot_1.jpeg.pdf"))
    unlink("plot_1.jpeg.pdf")
    expect_false(file.exists("plot_1.jpeg.pdf"))

    ## check that pdf suffix isn't added if not needed
    expect_is(diagnostic(out, K=1, pty="all", fname="plot_1.pdf"), "integer")
    expect_true(file.exists("plot_1.pdf"))
    unlink("plot_1.pdf")
    expect_false(file.exists("plot_1.pdf"))

    ## check that any suffix can be used when only 1 plot generated
    expect_is(diagnostic(out, K=1, pty="qq", fname="plot_qq.png"), "NULL")
    expect_true(file.exists("plot_qq.png"))
    unlink("plot_qq.png")
    expect_false(file.exists("plot_qq.png"))    

    ## check that error returned if invalid class fname specfied
    err_fname <- "fname must be NULL or a string"
    expect_error(diagnostic(out, K=1, pty="qq", fname=5000), err_fname)
})


test_that("plot accpets K parameter input", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, n_sim = 10, n_min = 30, alpha = 1)

    ## check that specifying vector K saves to pdf
    expect_is(diagnostic(out, K=1:3, pty='qq', fname="plot_123"), "integer")
    expect_true(file.exists("plot_123.pdf"))
    unlink("plot_123.pdf")
    expect_false(file.exists("plot_123.pdf"))

    ## check that ggplots are returned when only 1 plot specified
    plt_qq1 <- diagnostic(out, K=1, pty='qq')
    plt_qq2 <- diagnostic(out, K=2, pty='qq')
    expect_is(plt_qq1, "ggplot")
    expect_is(plt_qq2, "ggplot")

    ## check that K=2 and K=1 return different plots
    expect_gt(nrow(plt_qq1$data), nrow(plt_qq2$data))

    ## check that error returned with invalid K
    err_k <- "K must be a set of indices between 1 and n-1"
    expect_error(diagnostic(out, K=-100, pty='qq'), err_k)
    expect_error(diagnostic(out, K=1000, pty='qq'), err_k)    
})


test_that("plot accepts ci_idx parameter input", {
    ## simple dataset
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    out <- shc(dm, n_sim = 10, n_min = 30, alpha = 1,
               ci = c("2CI", "2CI"), null_alg = c("hclust", "2means"))

    ## check that ggplots are returned when only 1 plot specified
    plt_ci1 <- diagnostic(out, K=1, ci_idx = 1, pty='pvalue')
    plt_ci2 <- diagnostic(out, K=1, ci_idx = 2, pty='pvalue')
    expect_is(plt_ci1, "ggplot")
    expect_is(plt_ci2, "ggplot")

    ## check that points in plot match simulated CI values
    expect_equal(plt_ci1$data$x, out$ci_sim[1, , 1])
    expect_equal(plt_ci2$data$x, out$ci_sim[1, , 2])

    ## check that error returned if invalid ci_idx specified
    err_idx <- "invalid choice for ci_idx; ci_idx must be <="
    expect_error(diagnostic(out, K=1, ci_idx = 3, pty='pvalue'), err_idx)
})
