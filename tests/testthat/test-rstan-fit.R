test_that("rstan fit", {
  require(rstan)
  rstan::rstan_options("auto_write" = TRUE)
  m <- NULL
  stan_file_loc <- here::here("inst", "stan", "normal_group_means.stan")
  if (file.exists(stan_file_loc)) {
    m <- rstan::stan_model(file = stan_file_loc, save_dso = TRUE)
  }
  if (is.null(m)) {
    m <- rstan::stan_model(file = system.file(
      "stan", "normal_group_means.stan", package = "sbcrs"))
  }

  new_sbc_for_testing <- function(.n_obs, .n_groups, .n_types) {
    SBC$new(
      data = function(seed) {
        set.seed(seed + 10)
        n_obs <- .n_obs
        n_groups <- .n_groups
        n_types <- .n_types
        group <- sample.int(n_groups, size = n_obs, replace = TRUE)
        type <- sample.int(n_types, size = n_obs, replace = TRUE)
        data <- list(n_obs = n_obs, n_groups = n_groups, 
                     n_types = n_types, group = group, type = type)
        data
      },
      params = function(seed, data) {
        set.seed(seed + 20)
        sigma <- rexp(data$n_types)
        mu <- matrix(rnorm(data$n_groups * data$n_types, 0, 1),
          nrow = data$n_groups
        )
        nu <- rnorm(1)
        params <- list(sigma = sigma, mu = mu, nu = nu)
        params
      },
      modeled_data = function(seed, data, params) {
        set.seed(seed + 30)
        y_mean <- purrr::map2_dbl(data$group, data$type, ~ params$mu[.x, .y])
        y_sd <- params$sigma[data$type]
        modeled_data <- list(y = rnorm(data$n_obs, y_mean, y_sd))
        modeled_data
      },
      sampling = function(seed, data, params, modeled_data, iters) {
        samples <- rstan::sampling(m,
          data = c(data, modeled_data), seed = seed,
          chains = 1, iter = 2 * iters, warmup = iters
        )
        samples
      }
    )
  }

  check_param_len <- function(calib) {
    param_names <- names(calib$ranks)
    purrr::walk(param_names, ~ testthat::expect_equal(
      length(calib$ranks[[.x]]),
      length(calib$params[[.x]])
    ))
  }
  tests_for_sbc <- function(sbc, N, L, keep_stan_fit = TRUE) {
    testthat::expect_length(sbc$calibrations, 0)
    sbc$calibrate(N, L, keep_stan_fit)
    testthat::expect_length(sbc$calibrations, N)
    purrr::walk(sbc$calibrations, ~ check_param_len(.x))
    testthat::expect_s3_class(sbc$summary(), "data.frame")
    if (length(sbc$calibrations[[1]]$params$mu) > 0L)
      testthat::expect_s3_class(sbc$summary("mu"), "data.frame")
    gg <- sbc$plot()
    testthat::expect_s3_class(gg, "ggplot")
    if (length(sbc$calibrations[[1]]$params$mu) > 0L) {
      gg <- sbc$plot("mu")
      testthat::expect_s3_class(gg, "ggplot")
    }
    n_r <- sbc$calibrations[[1]]$ranks %>% unlist() %>% length()
    r <- sbc$ranks() %>% sort()
    testthat::expect_length(r, n_r * N)
    r0 <- purrr::map(sbc$calibrations, 'ranks') %>%
      purrr::map(unlist) %>%
      unlist() %>%
      unname() %>%
      sort()
    testthat::expect_equal(r, r0)
    if (keep_stan_fit) {
      testthat::expect_s4_class(sbc$calibrations[[1]]$samples, "stanfit")
    } else {
      testthat::expect_null(sbc$calibrations[[1]]$samples)
    }
    invisible(NULL)
  }
  sbc <- new_sbc_for_testing(0, 0, 0)
  tests_for_sbc(sbc, 4, 10)
  
  sbc <- new_sbc_for_testing(0, 0, 1)
  tests_for_sbc(sbc, 4, 10)
  
  sbc <- new_sbc_for_testing(0, 1, 1)
  tests_for_sbc(sbc, 4, 10)

  sbc <- new_sbc_for_testing(0, 3, 4)
  tests_for_sbc(sbc, N = 4, L = 10, keep_stan_fit = FALSE)

  sbc <- new_sbc_for_testing(0, 3, 4)
  tests_for_sbc(sbc, N = 4, L = 10, keep_stan_fit = FALSE)
  
  doParallel::registerDoParallel(cores = 2)
  sbc <- new_sbc_for_testing(0, 3, 4)
  tests_for_sbc(sbc, N = 4, L = 10, keep_stan_fit = FALSE)
  
})
