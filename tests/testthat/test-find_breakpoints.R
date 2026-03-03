test_that("find_breakpoints() works", {
  data(rbps)
  ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
  cemean <- find_breakpoints(ps, method = "cemean", Nmax = 5)
  ewma <- find_breakpoints(ps, method = "plateau", type = "ewma", lambda = 0.2)
  cusum <- find_breakpoints(ps, method = "plateau", type = "cusum")
  window <- find_breakpoints(ps, method = "window", window = 5)

  expect_true(inherits(cemean, "breakpoints"))

  expect_equal(cemean$end[1], 152)
  expect_equal(ewma$end[1], 157)
  expect_equal(cusum$end[1], 157)
  expect_true(window$end[1] >= 175)
})
