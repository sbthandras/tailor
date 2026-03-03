test_that("schoco works()", {
  data(rbps)
  data(adapters)
  amat <- adapter_matrix(adapters, ids = rbps$Core_ORF)
  out <- schoco(amat, k = 2)
  expect_equal(dim(out), c(2, 6))

  expect_error(schoco(adapters, 2))
  msg <- capture_error(schoco(adapters, 2))
  expect_equal(msg$message, paste0(
    "mat must be a matrix of class 'adapter_matrix'. ",
    "Use adapter_matrix() to create a compatible matrix."
  ))
})
