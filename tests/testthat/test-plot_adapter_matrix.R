test_that("plot() works with adapter matrices", {
  skip_on_ci()
  data(rbps)
  data(adapters)
  amat <- adapter_matrix(adapters)
  g1 <- plot(amat)
  clusters <- cluster_adapters(amat, k_min = 2, k_max = 2)
  g2 <- plot(amat, clusters = clusters)

  vdiffr::expect_doppelganger("amat-without-clusters", g1)
  vdiffr::expect_doppelganger("amat-with-clusters", g2)
})

test_that("plot() works with adapter matrices - covr hack", {
  data(rbps)
  data(adapters)
  amat <- adapter_matrix(adapters)
  g1 <- plot(amat)
  clusters <- cluster_adapters(amat, k_min = 2, k_max = 2)
  g2 <- plot(amat, clusters = clusters)

  expect_true(inherits(g1, "ggplot"))
  expect_true(inherits(g2, "ggplot"))
})
