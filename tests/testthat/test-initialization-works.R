test_that("initialization works", {
  sbc <- SBC$new(
    data = function(seed) {
      list(n = 10)
    },
    params = function(seed, data) {
      set.seed(seed + 10)
      list(mu = rnorm(data$n, 0, 1))
    },
    modeled_variable = function(seed, data, params) {
      set.seed(seed + 20)
      list(y = rnorm(data$n, params$mu, 1))
    },
    sampling = function(...) NULL
  )
  testthat::expect_setequal(class(sbc), c("SBC", "R6"))
  testthat::expect_type(sbc$calibrations, "list")
  testthat::expect_null(sbc$.data_fun)
  testthat::expect_length(sbc$calibrations, 0)
  testthat::expect_null(sbc$summary())
  testthat::expect_null(sbc$plot())
})
