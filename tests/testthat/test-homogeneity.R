test_that("homogeneity() works", {
  data(rbps)
  data(adapters)
  amat <- adapter_matrix(adapters, ids = rbps$Core_ORF)
  clusters <- cluster_adapters(amat, k_min = 2, k_max = 2)
  ids <- clusters |> dplyr::filter(cluster == "ACL 1") |> dplyr::pull(id)
  index <- which(rownames(amat) %in% ids)
  out <- homogeneity(mat = amat[index, index])

  expect_equal(out, 1)
})
