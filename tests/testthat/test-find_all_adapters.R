test_that("find_all_adapters() works", {
  data(rbps)
  adapters <- find_all_adapters(rbps$Core_ORF[1:3], data = rbps)

  expect_true(inherits(adapters, "adapter"))
  expect_equal(nrow(adapters), 3)
  expect_equal(adapters$pident, c(0.901, 0.940, 0.914))
})

test_that("find_all_adapters() work with multiple cores", {
  skip_on_cran()
  data(rbps)
  adapters <- find_all_adapters(rbps$Core_ORF[1:3], data = rbps, cores = 2)

  expect_true(inherits(adapters, "adapter"))
  expect_equal(nrow(adapters), 3)
  expect_equal(adapters$pident, c(0.901, 0.940, 0.914))
})
