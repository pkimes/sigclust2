context("shc plotting")


test_that("plot default call works properly", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, n_sim=10, n_min=30)

    ## check that ggplot is produced
    expect_is(plot(out), "ggplot")
})


test_that("plot accepts group labels", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, n_sim=10, n_min=30)

    ## different label classes
    grps_num <- c(rep(1, 50), rep(2, 50))
    grps_str <- c(rep("A", 50), rep("B", 50))
    grps_fac <- factor(c(rep("a", 50), rep("b", 50)))
    grps_bad <- c(rep("A", 100), rep("B", 20))

    plt_base <- plot(out)
    plt_num <- plot(out, groups = grps_num) 
    plt_str <- plot(out, groups = grps_str) 
    plt_fac <- plot(out, groups = grps_fac) 

    ## check that ggplots are generated for all options
    expect_is(plt_num, "ggplot")
    expect_is(plt_str, "ggplot")
    expect_is(plt_fac, "ggplot")

    ## check that label tiles are generated w/ groups specified
    layers_str <- sapply(plt_str$layers, function(x) class(x$geom)[1])
    expect_true("GeomTile" %in% layers_str)

    ## check that label tiles not generated in default plot
    layers_base <- sapply(plt_base$layers, function(x) class(x$geom)[1])
    expect_false("GeomTile" %in% layers_base)

    ## check that tile data labels match specified labels
    tile_idx <- which(layers_str == "GeomTile")[1]
    tile_dat <- plt_str$layers[[tile_idx]]$data
    tile_labs <- tile_dat$clusters[order(as.numeric(as.character(tile_dat$label)))]
    expect_equal(tile_labs, grps_str)
})


test_that("plot accepts use_labs parameter value", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, n_sim=10, n_min=30)

    ## note - assumes default, alpha=.05, FWER=TRUE
    plt_wlabs <- plot(out, use_labs = TRUE)
    plt_nolabs <- plot(out, use_labs = FALSE)

    ## check that ggplots are generated for all options
    expect_is(plt_wlabs, "ggplot")
    expect_is(plt_nolabs, "ggplot")

    ## check that only text labels generated w/ labs TRUE
    layers_wlabs <- sapply(plt_wlabs$layers, function(x) class(x$geom)[1])
    expect_true("GeomText" %in% layers_wlabs)
    expect_false("GeomTile" %in% layers_wlabs)

    ## check that text and tile labels not generated w/ labs FALSE
    layers_nolabs <- sapply(plt_nolabs$layers, function(x) class(x$geom)[1])
    expect_false("GeomText" %in% layers_nolabs)
    expect_false("GeomTile" %in% layers_wlabs)

    ## check that text data labels are just 1:n indices
    text_idx <- which(layers_wlabs == "GeomText")[1]
    text_dat <- plt_wlabs$layers[[text_idx]]$data
    text_labs <- sort(as.numeric(as.character(text_dat$label)))
    expect_equal(text_labs, 1:100)
})


test_that("plot accepts alpha, fwer parameter values", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    
    out_100 <- shc(dm, alpha = 1.00)
    out_010 <- shc(dm, alpha = 0.10)
    out_001 <- shc(dm, alpha = 0.01)

    plt_100_none <- plot(out_100, fwer = FALSE, alpha = 0.05)
    plt_100_001 <- plot(out_100, fwer = TRUE, alpha = 0.01)
    plt_010_001 <- plot(out_010, fwer = TRUE, alpha = 0.01)
    plt_001_none <- plot(out_001)
    
    ## check that ggplots are generated
    expect_is(plt_100_none, "ggplot")
    expect_is(plt_100_001, "ggplot")
    expect_is(plt_010_001, "ggplot")
    expect_is(plt_001_none, "ggplot")

    ## check that plots same when same 'fwer', 'alpha' specified to shc or plot
    gg_parts <- c("layers", "data", "scales", "mapping", "theme",
                  "coordinates", "facet", "labels")
    for (ipart in gg_parts) {
        expect_equal(plt_100_001[[ipart]], plt_001_none[[ipart]]) 
        expect_equal(plt_010_001[[ipart]], plt_001_none[[ipart]]) 
    }
    
    ## check that warning returned if plot 'alpha' smaller than shc 'alpha'
    warn_alpha <- "shc constructed using smaller alpha"
    expect_warning(plot(out_001, fwer = TRUE, alpha = 0.05), warn_alpha)

    ## check that warning refurned if plot 'fwer' disagrees w/ shc 'fwer'
    warn_fwer <- "shc constructed using alpha < 1, using fwer = TRUE"
    expect_warning(plot(out_001, fwer = FALSE, alpha = 0.05), warn_fwer)

    ## check that error returned if alpha outside of [0, 1]
    alpha_range_error <- "invalid choice for alpha; alpha must be 0 < alpha < 1"
    expect_error(plot(out_100, alpha = -.5), alpha_range_error)
    expect_error(plot(out_100, alpha = 1.5), alpha_range_error)
})


test_that("plot accepts ci_idx, ci_emp parameter values", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)

    out_100 <- shc(dm, alpha = 1, ci = c("2CI", "2CI"), null_alg = c("hclust", "2means"))
    out_010 <- shc(dm, alpha = 0.10, ci = c("2CI", "2CI"), null_alg = c("hclust", "2means"),
                   ci_idx = 2, ci_emp = TRUE)

    plt_100_ci1_t <- plot(out_100, ci_idx=1, ci_emp=TRUE)
    plt_100_ci1_f <- plot(out_100, ci_idx=1, ci_emp=FALSE)
    plt_100_ci2_t <- plot(out_100, ci_idx=2, ci_emp=TRUE)
    plt_100_ci2_f <- plot(out_100, ci_idx=2, ci_emp=FALSE)

    plt_100_010_ci2_t <- plot(out_100, alpha=0.10, ci_idx=2, ci_emp=TRUE)
    plt_010_ci2_t <- plot(out_010, alpha=0.10)

    ## check that ggplots are generated
    expect_is(plt_100_ci1_t, "ggplot")
    expect_is(plt_100_ci1_f, "ggplot")
    expect_is(plt_100_ci2_t, "ggplot")
    expect_is(plt_100_ci2_f, "ggplot")
    expect_is(plt_010_ci2_t, "ggplot")

    ## check that plots same when same 'fwer', 'alpha' specified to shc or plot
    gg_parts <- c("layers", "data", "scales", "mapping", "theme",
                  "coordinates", "facet", "labels")
    for (ipart in gg_parts) {
        expect_equal(plt_100_010_ci2_t[[ipart]],
                     plt_010_ci2_t[[ipart]]) 
    }

    ## check that error returned when invalid ci_idx specified
    err_idx <- "invalid choice for ci_idx"
    expect_error(plot(out_010, ci_idx = 3, ci_emp = TRUE), err_idx)
    expect_error(plot(out_010, ci_idx = 3, ci_emp = FALSE), err_idx)

    ## check that warning returned when different ci_idx specified
    warn_idx <- "shc constructed using alpha < 1"
    expect_warning(plot(out_010, ci_idx = 2, ci_emp = FALSE), warn_idx)

    ## check that warning returned when different ci_emp specified
    warn_emp <- "shc constructed using alpha < 1"
    expect_warning(plot(out_010, ci_idx = 1, ci_emp = TRUE), warn_emp)

})


test_that("plot accpets hang parameter values", {
    ## simple dataset
    set.seed(100)
    dm <- matrix(rnorm(200), ncol=2, nrow=100)
    out <- shc(dm, alpha = 1)

    ## check that specifying 'hang' doesn't cause error
    expect_is(plot(out, hang = -1), "ggplot")
    expect_is(plot(out, hang = 0), "ggplot")
    expect_is(plot(out, hang = 0.5), "ggplot")
})
