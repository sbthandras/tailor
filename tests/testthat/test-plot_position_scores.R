test_that("plot() works with position scores", {
  skip_on_ci()
  data(rbps)
  ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
  g <- plot(ps)

  vdiffr::expect_doppelganger("position-scores", g)
})
