test_that("adapter_matrix() works", {
  data(adapters)
  amat <- adapters |> adapter_matrix()

  expect_true(inherits(amat, "adapter_matrix"))
  expect_equal(dim(amat), c(20,20))
  expect_equal(unique(diag(amat)), 1)

  index_x <- which(rownames(amat) == "MN395291-1")
  index_y <- which(colnames(amat) == "ON513429-1")
  expect_equal(amat[index_x, index_y], 0.876)
})

test_that("adapter_matrix() work with multiple cores", {
  skip_on_cran()
  data(adapters)
  amat <- adapters |> adapter_matrix(cores = 2)

  expect_true(inherits(amat, "adapter_matrix"))
  expect_equal(dim(amat), c(20,20))
  expect_equal(unique(diag(amat)), 1)

  index_x <- which(rownames(amat) == "MN395291-1")
  index_y <- which(colnames(amat) == "ON513429-1")
  expect_equal(amat[index_x, index_y], 0.876)
})

test_that("adapter_matrix() fails when input is not an adapter df", {
  data(rbps)

  expect_error(rbps |> adapter_matrix())
  msg <- capture_error(rbps |> adapter_matrix())

  expect_equal(msg$message, paste0(
    "'adapters' must be a data frame of class 'adapter'. ",
    "Use find_adapter() to create a compatible data frame."
  ))
})

test_that("adapter_matrix() fails for id mismatch", {
  data(adapters)
  expect_error(adapters |> adapter_matrix(ids = "kutya"))

  msg <- capture_error(adapters |> adapter_matrix(ids = "kutya"))
  expect_equal(
    msg$message,
    "The following ids were not found in the 'adapter' data frame: kutya"
  )

  ids <- c(adapters$pattern_id, adapters$subject_id) |> unique()
  ids <- ids[-1]
  expect_error(adapters |> adapter_matrix(ids = ids))

  msg <- capture_error(adapters |> adapter_matrix(ids = ids))
  expect_equal(
    msg$message,
    "The following ids were not found in the 'ids' vector: MN395291-1"
  )

})


