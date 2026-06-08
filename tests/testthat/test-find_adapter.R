test_that("find_adapter() works", {
  data(rbps)
  ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
  bps <- find_breakpoints(ps, Nmax = 5)
  ada <- find_adapter(bps)

  expect_true(inherits(ada, "adapter"))
  expect_equal(nrow(ada), 1)
  expect_equal(ada$pident, 0.901)
})

test_that("find_adapter() fails when input is not a breakpoints object", {
  data(rbps)
  ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
  expect_error(find_adapter(ps))
  msg <- capture_error(find_adapter(ps))
  expect_equal(msg$message, paste0(
    "bps must be a data frame of class 'breakpoints'. ",
    "Use find_breakpoints() to create a compatible data frame."
  ))
})

test_that("find_adapter() returns NAs if no adapter found", {
  data(rbps)
  ps <- position_scores("MW056503-1", "OX335376-1", data = rbps)
  bps <- find_breakpoints(ps, Nmax = 5)
  ada <- find_adapter(bps)

  expect_true(inherits(ada, "adapter"))
  expect_equal(nrow(ada), 1)
  expect_equal(ada$start, NA_integer_)
  expect_equal(ada$end, NA_integer_)
  expect_equal(ada$mean_score, NA_real_)
  expect_equal(ada$pident, NA_real_)
})
