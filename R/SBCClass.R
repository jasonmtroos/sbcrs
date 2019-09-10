


#' SBC class
#'
#' @description Generator class for creating new instances of the `SBC` R6 class.
#'
#' @format An R6 object of type `SBC`
#'
#' @details See [`vignette('intro-to-sbc', package = 'sbcrs')`](../doc/intro-to-sbc.html) for
#' an accessible introduction to simulation-based calibration.
#'
#' @section Methods:
#' \describe{
#' \item{`$new(data, params, modeled_data, sampling)`}{Create a new SBC object, passing in functions to generate data and parameters, and draw samples.
#'     \describe{
#'         \item{\code{data = function(seed) {}}}{A function with signature `function(seed)` that returns a named list.}
#'         \item{\code{params = function(seed, data) {}}}{A function with signature `function(seed, data)` that returns a named list.}
#'         \item{\code{modeled_data = function(seed, data, params) {}}}{A function with signature `function(seed, data, params)` that returns a named list.}
#'         \item{\code{sampling = function(seed, data, params, modeled_data, iters) {}}}{A function with signature `function(seed, data, params, modeled_data, iters)` that returns a `stanfit` object run for `iters` sampling iterations.}
#'     }}
#' \item{`$calibrate(N, L, keep_stan_fit = TRUE)`}{
#'     Run the calibration procedure.
#'     \describe{
#'         \item{`N`}{The number of times to simulate parameters and recover via MCMC.}
#'         \item{`L`}{The number of MCMC samples to retain when calculating rank statistics.}
#'         \item{`keep_stan_fit = TRUE`}{If `TRUE` (the default), then the `stan_fit` objects returned by the `sampling` function are retained.}
#'     }
#' }
#' \item{`$summary(var = NULL)`}{Summarize results of a previous calibration. Optionally specify a parameter `var`.}
#' \item{`$plot(var = NULL)`}{Plot a histogram of ranks from a previous calibration. Optionally specify a parameter `var`.}
#' }

#'
#' @section Fields:
#' \describe{
#' \item{`$calibrations`}{A list of `N` calibrations created by calling `$calibrate`. Each item is list with the following named elements:
#'     \describe{
#'        \item{data}{Output from call to the `data` function.}
#'        \item{params}{Output from call to the `params` function.}
#'        \item{modeled_data}{Output from call to the `modeled_data` function.}
#'        \item{samples}{A stanfit object returned from call to the `sampling` function if `keep_stan_fit = TRUE`, otherwise NULL.}
#'        \item{ranks}{A named list matching items in `params`. Values express the number of samples out of a maximum `L` with `sampled$var < param$var`, where `sampled$var` indicates a vector of `L` samples of parameter `var`.}
#'        \item{n_eff}{The smallest effective sample size for any parameter in any of the `stan_fit` objects.}
#'        \item{iters}{The number of MCMC iterations (before thinning) in the `stan_fit` object.}
#'        \item{seed}{The value of `seed` passed to the data, parameter, and sampling functions.}
#'     }
#' }
#' }
#'
#'
#'
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' sbc <- SBC$new(
#'   data = function(seed) {
#'     list(n = 10)
#'   },
#'   params = function(seed, data) {
#'     list(mu = rnorm(1))
#'   },
#'   modeled_data = function(seed, data, params) {
#'     list(y = rnorm(data$n, mu, 1))
#'   },
#'   sampling = function(seed, data, params, modeled_data) {
#'     stan_object <- NULL # usually a call to rstan::sampling()
#'     stan_object
#'   }
#' )
#' }
#' @export
SBC <- R6::R6Class(
  classname = "SBC",
  private = list(
    .data_fun = NULL,
    .params_fun = NULL,
    .modeled_data_fun = NULL,
    .sampling_fun = NULL,
    .N = NULL,
    .L = NULL,
    .new_calibration =
      function(seed, L, keep_stan_fit) {
        data <- private$.data_fun(seed)
        params <- private$.params_fun(seed, data)
        modeled_data <- private$.modeled_data_fun(seed, data, params)

        iters <- 1000
        n_eff <- as.double(L) - .0001

        while (n_eff < L) {
          iters <- round(iters * L / n_eff)

          samples <- private$.sampling_fun(
            seed, data, params, modeled_data, iters)
          # from rstan:::print.stan_fit:
          n_eff <- min(round(rstan::summary(
            samples)$summary[, "n_eff"], 0), na.rm = TRUE) 

          if (length(samples@sim) == 0) {
            stop(paste0("Sampling problem"))
          }
          if (n_eff < L) {
            warning("n_eff is smaller than L. Resampling")
          }
        }

        # from rstan:::print.stan_fit:
        n_samples <- samples@sim$n_save - samples@sim$warmup2 
        thin <- round(seq.int(1, n_samples, length.out = L))


        ranks_for_param <- function(par, par_name) {
          s <- rstan::extract(samples, par_name)
          if (length(s) == 0L) {
            return(NULL)
          }
          s <- s[[1]]

          dim_p <- dim(par) %||% length(par)
          if (sum(dim_p) == 0L) {
            stop("Zero-length parameters are not allowed")
          }

          if (length(dim_p) == 1L && dim_p[1] == 0L) {
            # scalar parameter
            q <- sum(s[thin] < par)
          } else if (length(dim_p) == 1L) {
            # vector parameter
            q <- colSums(matrix(s, ncol = dim_p)[thin, ] < 
                           (rep(1, L) %*% t(par)))
          } else {
            # matrix++ parameter
            q <- plyr::aaply(s[thin, , ], .margins = 1, 
                             .fun = function(x) as.vector(x < par))
            dim_q <- c(length(thin), dim_p)
            dim(q) <- dim_q
            q <- unname(colSums(q))
          }
        }
        ranks <- purrr::imap(params, ~ ranks_for_param(.x, .y))

        if (!keep_stan_fit) {
          samples <- NULL
        }
        list(
          data = data,
          params = params,
          modeled_data = modeled_data,
          samples = samples,
          ranks = ranks,
          n_eff = n_eff,
          iters = iters,
          seed = seed
        )
      },
    .summary_data = function(var = NULL) {
      q <- purrr::map(self$calibrations, "ranks")
      q <- purrr::map(q, ~.x[!purrr::map_lgl(.x, is.null)])
      q <- purrr::imap(
        purrr::transpose(q),
        ~ matrix(purrr::simplify(.x),
          nrow = length(q),
          byrow = TRUE
        )
      )

      if (is.null(var)) {
        qd <- tibble::tibble(r = as.vector(unlist(q)))
        N <- nrow(qd)
      } else {
        qd <- tidyr::gather(dplyr::mutate(dplyr::tbl_df(as.data.frame(q[var])),
          sim = dplyr::row_number()
        ), var, r, -sim)
        N <- private$.N
      }
      list(N = N, qd = qd)
    }
  ),
  public = list(
    initialize =
      function(data = function(seed) list(),
                     params = function(seed, data) list(),
                     modeled_data = function(seed, data, params) list(),
                     sampling) {
        private$.data_fun <- data
        private$.params_fun <- params
        private$.modeled_data_fun <- modeled_data
        private$.sampling_fun <- sampling
        invisible(self)
      },
    calibrations = list(),
    calibrate = function(N, L, keep_stan_fit = TRUE) {
      stopifnot(N > 0L && L > 1L)
      if (foreach::getDoParWorkers() == 1) {
        self$calibrations <- purrr::map(
          seq_len(N), ~ private$.new_calibration(
            seed = .x, L = L, keep_stan_fit = keep_stan_fit))
      } else {
        `%dopar%` <- foreach::`%dopar%`
        self$calibrations <- foreach::foreach(
          seed = seq_len(N)) %dopar%
          private$.new_calibration(
            seed = seed, L = L, keep_stan_fit = keep_stan_fit)
      }

      smallest_n_eff <- min(purrr::map_dbl(
        self$calibrations, "n_eff"))
      if (L > smallest_n_eff) {
        warning("L is too large. Smallest n_eff = ", round(n_eff))
      }

      private$.N <- N
      private$.L <- L

      invisible(self)
    },
    summary = function(var = NULL) {
      if (is.null(private$.N)) {
        return(invisible(NULL))
      }

      sd <- private$.summary_data(var)
      N <- sd$N
      qd <- sd$qd
      L <- private$.L
      iqs <-
        purrr::map(c(.5, 1 - 1 / L), ~ binom_lims(.x, N, L)) %>%
        dplyr::bind_rows()

      if (is.null(var)) {
        dd <-
          qd %>%
          dplyr::group_by(r) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::mutate(var = "")
      } else {
        dd <-
          qd %>%
          dplyr::group_by(r, var) %>%
          dplyr::summarise(n = dplyr::n())
      }

      dd2 <-
        dd %>%
        tidyr::crossing(dplyr::distinct(iqs)) %>%
        dplyr::mutate(inside = n >= lo & n <= hi, 
                      outside = n < lo | n > hi) %>%
        dplyr::group_by(var, iq) %>%
        dplyr::summarise(inside = sum(inside), 
                         outside = sum(outside)) %>%
        dplyr::mutate(inside = inside / (inside + outside), 
                      outside = 1 - inside) %>%
        tidyr::gather(key, actual, -iq, -var) %>%
        dplyr::mutate(expected = dplyr::case_when(
          key == "inside" ~ iq,
          TRUE ~ 1 - iq
        )) %>%
        dplyr::arrange(key, iq, var) %>%
        dplyr::filter(key == "outside") %>%
        dplyr::select(var, iq, expected.outside = expected, 
                      actual.outside = actual) %>%
        dplyr::ungroup()

      purrr::walk(
        split(dd2, dd2$var),
        function(x) {
          cat(paste0("
", x$var[1], "
"))
          print(format.data.frame(dplyr::select(x, -var)), 
                row.names = FALSE)
        }
      )
      invisible(dd2)
    },
    plot = function(var = NULL) {
      if (is.null(private$.N)) {
        return(invisible(NULL))
      }

      L <- private$.L
      sd <- private$.summary_data(var)
      N <- sd$N
      qd <- sd$qd
      if (is.null(var)) {
        facets <- list(ggplot2::facet_null())
      } else {
        facets <- list(ggplot2::facet_wrap(~var))
      }

      iqs <-
        c(.5, 1 - 1 / L) %>%
        purrr::map(~ binom_lims(.x, N, L)) %>%
        dplyr::bind_rows() %>%
        tidyr::crossing(r = c(0, L))

      plotfun <- function() {
        ggplot2::ggplot(iqs, aes(
          x = r, ymin = lo, ymax = hi, alpha = factor(iq))) +
          ggplot2::geom_ribbon(fill = "black") +
          ggplot2::stat_count(
            ggplot2::aes(ymin = NULL, ymax = NULL),
            data = qd, alpha = 1, size = 3,
            fill = "white", colour = "black", 
            geom = "point", shape = 21
          ) +

          ggplot2::scale_alpha_discrete(range = c(.2, .1)) +
          ggplot2::theme_minimal() +
          facets +
          ggplot2::labs(alpha = "IQR", x = stringr::str_glue(
            "Quantile rank ({L} MCMC draws)"), y = "Realizations")
      }
      qpf <- purrr::quietly(plotfun)
      qpf <- qpf()$result
      qpf
    }
  ),
  inherit = NULL,
  lock_objects = TRUE,
  class = TRUE,
  portable = TRUE,
  lock_class = FALSE,
  cloneable = TRUE
)
