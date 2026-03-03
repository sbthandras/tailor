test_that("completeness() works", {
  data(rbps)
  data(adapters)
  amat <- adapter_matrix(adapters)
  clusters <- cluster_adapters(amat, k_min = 2, k_max = 2)
  ids <- clusters |> dplyr::filter(cluster == "ACL 1") |> dplyr::pull(id)
  index <- which(rownames(amat) %in% ids)
  out <- completeness(mat = amat, index = index)

  expect_equal(out, 1)
})
