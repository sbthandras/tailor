test_that("cluster_adapters() works", {
  data(rbps)
  data(adapters)
  amat <- adapter_matrix(adapters)

  clusters1 <- cluster_adapters(amat, k_min = 1)
  expect_true(inherits(clusters1, "data.frame"))
  expect_equal(dim(clusters1), c(20,2))
  expect_equal(clusters1 |> dplyr::pull(cluster) |> unique() |> length(), 2)

  clusters2 <- cluster_adapters(amat, k_min = 2, k_max = 5)
  expect_true(inherits(clusters2, "data.frame"))
  expect_equal(dim(clusters2), c(20,2))
  expect_equal(clusters2 |> dplyr::pull(cluster) |> unique() |> length(), 2)
})

test_that("cluster_adapters() works with multiple cores", {
  skip_on_cran()
  data(rbps)
  data(adapters)
  amat <- adapter_matrix(adapters)

  clusters1 <- cluster_adapters(amat, k_min = 1, cores = 2)
  expect_true(inherits(clusters1, "data.frame"))
  expect_equal(dim(clusters1), c(20,2))
  expect_equal(clusters1 |> dplyr::pull(cluster) |> unique() |> length(), 2)
})

test_that("cluster_adapters() fails when input is not an adapter matrix", {
  data(rbps)
  expect_error(cluster_adapters(rbps, k_min = 2, k_max = 5))

  msg <- capture_error(cluster_adapters(rbps, k_min = 2, k_max = 5))
  expect_equal(msg$message, paste0(
   "mat must be a matrix of class 'adapter_matrix'. ",
   "Use adapter_matrix() to create a compatible matrix."
  ))
})
