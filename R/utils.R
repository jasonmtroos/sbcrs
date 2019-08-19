binom_lims <-
  function(iq, N, L) {
    .x <- stats::qbinom(.5 + c(-iq/2, iq/2), N, 1/(L+1))
    tibble::tibble(lo = .x[1], hi = .x[2], iq = iq)
  }


`%||%` <- function(a, b) {
  a[is.null(a)] <- b[is.null(a)]
  a
}
