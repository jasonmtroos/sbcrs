
#' @export
SBC <- R6::R6Class(
  classname = 'SBC',
  private = list(
    .data_fun = NULL,
    .params_fun = NULL,
    .modeled_variable_fun = NULL,
    .sampling_fun = NULL,
    .N = NULL,
    .L = NULL,
    .new_calibration = 
      function(seed, L) {
        data <- private$.data_fun(seed)
        params <- private$.params_fun(seed, data)
        modeled_variable <- private$.modeled_variable_fun(seed, data, params)
        
        iters <- 1000
        n_eff <- as.double(L) - .0001
        
        while (n_eff < L) {
          iters <- round(iters * L / n_eff)
          
          samples <- private$.sampling_fun(seed, data, params, modeled_variable, iters)
          n_eff <- round(min((as.data.frame(rstan::summary(samples)$summary))$n_eff, na.rm = TRUE))
          
          if (length(samples@sim) == 0) {
            stop(paste0('Sampling problem'))
          }
          if (n_eff < L) {
            warning('n_eff is smaller than L. Resampling')
          }
        }
        
        n_samples <- samples@stan_args[[1]]$iter - samples@stan_args[[1]]$warmup
        thin <- round(seq.int(1, n_samples, length.out = L))
        
        
        q_for_p <- function(par, par_name) {
          s <- rstan::extract(samples, par_name)
          if (length(s) == 0L)
            return(NULL)
          s <- s[[1]]
          
          dim_p <- dim(par) %||% length(par)
          if (sum(dim_p) == 0L)
            stop('Zero-length parameters are not allowed')
          
          if (length(dim_p) == 1L && dim_p[1] == 0L) {
            # scalar parameter
            q <- sum(s[thin] < par)
          } else if (length(dim_p) == 1L) {
            # vector parameter
            q <- colSums(matrix(s, ncol = dim_p)[thin,] < (rep(1, L) %*% t(par)))
          } else {
            # matrix++ parameter
            q <- plyr::aaply(s[thin,,], .margins = 1, .fun = function(x) as.vector(x < par))
            if (is.null(dim(q))) {
              q <- array(q, dim = dim(s)[-1])
            }
            q <- unname(colSums(q))
          }
        }
        quantiles <- purrr::imap(params, ~q_for_p(.x, .y))
        
        list(data = data, 
             params = params, 
             modeled_variable = modeled_variable,
             samples = samples,
             quantiles = quantiles,
             n_eff = n_eff,
             iters = iters,
             seed = seed)
      },
    .clean_quantiles = function() {
      quantiles <- purrr::map(self$calibrations, 'quantiles')
      purrr::imap(purrr::transpose(quantiles), 
                  ~ matrix(purrr::simplify(.x), 
                           nrow = length(quantiles), 
                           byrow = TRUE))
    },
    .summary_data = function(var = NULL) {
      L <- private$.L
      q <- private$.clean_quantiles()
      if (is.null(var)) {
        qd <- tibble::tibble(r = as.vector(unlist(q)))
        N <- nrow(qd)
      } else {
        qd <- tidyr::gather(dplyr::mutate(dplyr::tbl_df(as.data.frame(q[var])),
                                          sim = dplyr::row_number()), var, r, -sim)
        N <- private$.N
      }
      list(N = N, qd = qd)
    }
  ),
  public = list(
    initialize = 
      function(data = NULL, 
               params = NULL, 
               modeled_variable = NULL,
               sampling = NULL) {
        private$.data_fun <- data
        private$.params_fun <- params
        private$.modeled_variable_fun <- modeled_variable
        private$.sampling_fun <- sampling
        invisible(self)
      },
    calibrations = list(),
    calibrate = function(N, L) {
      stopifnot(N > 0L && L > 1L)
      if (foreach::getDoParWorkers() == 1) {
        self$calibrations <- purrr::map(seq_len(N), ~private$.new_calibration(seed = .x, L = L))
      } else {
        `%dopar%` <- foreach::`%dopar%`
        self$calibrations <- foreach::foreach(seed = seq_len(N)) %dopar% 
          private$.new_calibration(seed = seed, L = L)
      }
      
      smallest_n_eff <- min(purrr::map_dbl(self$calibrations, 'n_eff'))
      if (L > smallest_n_eff)
        warning('L is too large. Smallest n_eff = ', round(n_eff))
      
      private$.N <- N
      private$.L <- L
      
      invisible(self)
    },
    summary = function(var = NULL) {
      if (is.null(private$.N))
        return(invisible(NULL))
      
      sd <- private$.summary_data(var)
      N <- sd$N
      qd <- sd$qd
      L <- private$.L
      iqs <-
        purrr::map(c(.5, 1-1/L), ~binom_lims(.x, N, L)) %>%
        dplyr::bind_rows()
      
      if (is.null(var)) {
        dd <-
          qd %>%
          dplyr::group_by(r) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::mutate(var = '')
      } else {
        dd <-
          qd %>%
          dplyr::group_by(r, var) %>%
          dplyr::summarise(n = dplyr::n())
      }
      
      dd2 <-
        dd %>%
        tidyr::crossing(dplyr::distinct(iqs)) %>%
        dplyr::mutate(inside = n >= lo & n <= hi, outside = n < lo | n > hi) %>%
        dplyr::group_by(var, iq) %>%
        dplyr::summarise(inside = sum(inside), outside = sum(outside)) %>%
        dplyr::mutate(inside = inside / (inside + outside), outside = 1 - inside) %>%
        tidyr::gather(key, actual, -iq, -var) %>%
        dplyr::mutate(expected = dplyr::case_when(key == 'inside' ~ iq,
                                                  TRUE ~ 1 - iq)) %>%
        dplyr::arrange(key, iq, var) %>%
        dplyr::filter(key == 'outside') %>%
        dplyr::select(var, iq, expected.outside = expected, actual.outside = actual) %>%
        dplyr::ungroup()
      
      purrr::walk(split(dd2, dd2$var), 
                  function(x) {
                    cat(paste0('\n', x$var[1], '\n'))
                    print(format.data.frame(dplyr::select(x, -var)), row.names = FALSE)
                  })
      invisible(dd2)
    },
    plot = function(var = NULL) {
      if (is.null(private$.N))
        return(invisible(NULL))
      
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
        c(.5, 1-1/L) %>%
        purrr::map(~binom_lims(.x, N, L)) %>%
        dplyr::bind_rows() %>%
        tidyr::crossing(r = c(0, L))
      
      plotfun <- function() {
        ggplot2::ggplot(iqs, aes(x = r, ymin = lo, ymax = hi, alpha = factor(iq))) +
          ggplot2::geom_ribbon(fill = 'black') +
          ggplot2::stat_count(ggplot2::aes(ymin = NULL, ymax = NULL), data = qd, alpha = 1, size = 3,
                              fill = 'white', colour = 'black', geom = 'point', shape = 21) +
          
          ggplot2::scale_alpha_discrete(range = c(.2, .1)) +
          ggplot2::theme_minimal() +
          facets +
          ggplot2::labs(alpha = 'IQR', x = stringr::str_glue('Quantile rank ({L} MCMC draws)'), y = 'Realizations')
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
  cloneable = TRUE)
             
             
             
