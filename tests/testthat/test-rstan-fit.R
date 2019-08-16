test_that("rstan fit", {
  library(rstan)
  
  rstan::rstan_options("auto_write" = TRUE)
  m <- NULL
  if (interactive()) {
    stan_file_loc <- here::here('inst', 'stan', 'normal_group_means.stan')
    if (file.exists(stan_file_loc)) {
      m <- stan_model(file = stan_file_loc, save_dso = TRUE)
    }
  }
  if (is.null(m)) {
    m <- stan_model(file = system.file('stan', 'normal_group_means.stan', package = 'sbcrs'))
  }
  
  new_sbc_for_testing <- function(.n_obs, .n_groups, .n_types) {
    stopifnot(.n_groups > 0L && .n_types > 0L)
    sbc <- SBC$new(
      data = function(seed) {
        set.seed(seed + 10)
        n_obs <- .n_obs
        n_groups <- .n_groups
        n_types <- .n_types
        group <- sample.int(n_groups, size = n_obs, replace = TRUE)
        type <- sample.int(n_types, size = n_obs, replace = TRUE)
        data <- list(n_obs = n_obs, n_groups = n_groups, n_types = n_types, group = group, type = type)
        data
      },
      params = function(seed, data) {
        set.seed(seed + 20)
        sigma <- rexp(data$n_types)
        mu <- matrix(rnorm(data$n_groups * data$n_types, 0, 1), 
                     nrow = data$n_groups)
        nu  <- rnorm(1)
        params <- list(sigma = sigma, mu = mu, nu = nu)
        params
      },
      modeled_variable = function(seed, data, params) {
        set.seed(seed + 30)
        y_mean <- purrr::map2_dbl(data$group, data$type, ~params$mu[.x, .y])
        y_sd <- params$sigma[data$type]
        modeled_variable <- list(y = rnorm(data$n_obs, y_mean, y_sd))
        modeled_variable
      },
      sampling = function(seed, data, params, modeled_variable, iters) {
        stan_sample_shhhh <- purrr::quietly(sampling)
        result <- stan_sample_shhhh(m, data = c(data, modeled_variable), seed = seed,
                                    chains = 1, iter = 2 * iters, warmup = iters)
        samples <- result$result
        samples
      })
  }
  
  check_param_len <- function(calib) {
    param_names <- names(calib$quantiles)
    purrr::walk(param_names, ~testthat::expect_equal(length(calib$quantiles[[.x]]), 
                                                     length(calib$params[[.x]])))
  }
  tests_for_sbc <- function(sbc, N, L) {
    testthat::expect_length(sbc$calibrations, 0)
    sbc$calibrate(N, L)
    testthat::expect_length(sbc$calibrations, N)
    purrr::walk(sbc$calibrations, ~check_param_len(.x))
    testthat::expect_s3_class(purrr:::quietly(sbc$summary)()$result, 'data.frame')
    gg <- sbc$plot()
    testthat::expect_s3_class(gg, 'ggplot')
    invisible(NULL)
  }
  sbc <- new_sbc_for_testing(0, 1, 1)
  tests_for_sbc(sbc, 1, 2)
  
  sbc <- new_sbc_for_testing(0, 3, 4)
  tests_for_sbc(sbc, N = 4, L = 10)
})
