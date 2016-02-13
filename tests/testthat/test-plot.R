library("sigclust2")
context("shc output plotting")

test_that("plot default call works properly", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, alpha = 1)
    
    expect_is(plot(out), "ggplot")
})


test_that("plot groups input is properly parsed", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, alpha = 1)

    grps_num <- c(rep(1, 50), rep(2, 50))
    grps_str <- c(rep("A", 50), rep("B", 50))
    grps_fac <- factor(c(rep("a", 50), rep("b", 50)))

    expect_is(plot(out, groups = grps_num), "ggplot")
    expect_is(plot(out, groups = grps_str), "ggplot")
    expect_is(plot(out, groups = grps_fac), "ggplot")
})


test_that("plot use_labs input is properly parsed", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, alpha = 1)

    expect_is(plot(out, use_labs = TRUE), "ggplot")
    expect_is(plot(out, use_labs = FALSE), "ggplot")
})


test_that("plot alpha, fwer inputs are properly parsed", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    
    fwer_false_warning <- "shc constructed using alpha < 1, using fwer = TRUE"
    alpha_small_warning <- "shc constructed using smaller alpha than specified to plot"
    alpha_range_error <- "invalid choice for alpha; alpha must be 0 < alpha < 1"    

    out1 <- shc(dm, alpha = 1.00)
    expect_is(plot(out1, fwer = TRUE, alpha = 0.05), "ggplot")
    expect_is(plot(out1, fwer = FALSE, alpha = 0.05), "ggplot")

    out2 <- shc(dm, alpha = 0.10)
    expect_is(plot(out2, fwer = TRUE, alpha = 0.05), "ggplot")
    expect_warning(plot(out2, fwer = FALSE, alpha = 0.05), fwer_false_warning)

    out3 <- shc(dm, alpha = 0.01)
    expect_warning(plot(out3, fwer = TRUE, alpha = 0.05), alpha_small_warning)
    expect_warning(plot(out3, fwer = FALSE, alpha = 0.05), fwer_false_warning)
    
    expect_error(shc(dm, alpha = -.5), alpha_range_error)
    expect_error(shc(dm, alpha = 1.5), alpha_range_error)    
})


test_that("plot hang input is properly parsed", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, alpha = 1)

    expect_is(plot(out, hang = -1), "ggplot")
    expect_is(plot(out, hang = 0), "ggplot")
    expect_is(plot(out, hang = 0.5), "ggplot")
})


test_that("plot ci_idx, ci_emp inputs are properly parsed", {
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out1 <- shc(dm, alpha = 1, ci = c("2CI", "2CI"), null_alg = c("hclust", "2means"))
    out2 <- shc(dm, alpha = 0.10, ci = c("2CI", "2CI"), null_alg = c("hclust", "2means"))

    ci_emp_warning <- "shc constructed using alpha < 1"#, using ci_emp from in_args(shc)"
    ci_idx_warning <- "shc constructed using alpha < 1"#, using ci_idx from in_args(shc)"
    ci_idx_error <- "invalid choice for ci_idx"#; ci_idx must be < length(ci)"
    
    expect_is(plot(out1, ci_idx = 1, ci_emp = TRUE), "ggplot")
    expect_is(plot(out1, ci_idx = 1, ci_emp = FALSE), "ggplot")
    expect_is(plot(out1, ci_idx = 2, ci_emp = TRUE), "ggplot")
    expect_is(plot(out1, ci_idx = 2, ci_emp = FALSE), "ggplot")

    expect_error(plot(out1, ci_idx = 3, ci_emp = TRUE), ci_idx_error)
    expect_error(plot(out1, ci_idx = 3, ci_emp = FALSE), ci_idx_error)

    expect_is(plot(out2), "ggplot")
    expect_is(plot(out2, ci_idx = 1, ci_emp = FALSE), "ggplot")
    expect_warning(plot(out2, ci_idx = 2, ci_emp = FALSE), ci_idx_warning)
    expect_warning(plot(out2, ci_idx = 1, ci_emp = TRUE), ci_emp_warning)
})
