test_that("rstan sample priors from stan", {
  library(rstan)
  rstan::rstan_options("auto_write" = TRUE)
  m <- NULL
  if (interactive()) {
    stan_file_loc <- here::here('inst', 'stan', 'rstan_sbc_example_modified.stan')
    if (file.exists(stan_file_loc)) {
      m <- stan_model(file = stan_file_loc, save_dso = TRUE)
    }
  }
  if (is.null(m)) {
    m <- stan_model(file = system.file('stan', 'normal_group_means.stan', package = 'sbcrs'))
  }
  
  sbc <- SBC$new(
      data = function(seed) {
        N <- 10
        a <- 2
        b <- 2
        data <- list(N = N, a = a, b = b)
        data
      },
      params = function(seed, data) {
        stan_sample_shhhhh <- purrr::quietly(sampling)
        pars <- extract(stan_sample_shhhhh(m, pars = 'pi_', chains = 1,
                                           seed = seed + 10,
                                   data = c(data, list(y = data$N)), iter = 1)$result)
        params <- list(pi = pars$pi_)
        params
      },
      modeled_variable = function(seed, data, params) {
        set.seed(seed + 20)
        modeled_variable <- list(y = rbinom(1, data$N, params$pi))
        modeled_variable
      },
      sampling = function(seed, data, params, modeled_variable, iters) {
        stan_sample_shhhh <- purrr::quietly(sampling)
        result <- stan_sample_shhhh(m, data = c(data, modeled_variable), seed = seed,
                                    chains = 1, iter = 2 * iters, warmup = iters,
                                    include = FALSE, pars = 'pi_')
        samples <- result$result
        samples
      })
  testthat::expect_length(sbc$calibrations, 0)
  sbc$calibrate(2, 10)
  testthat::expect_length(sbc$calibrations, 2)
})
