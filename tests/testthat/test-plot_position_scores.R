test_that("plot() works with position scores", {
  skip_on_ci()
  data(rbps)
  ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
  g1 <- plot(ps)
  g2 <- plot(ps, method = "cemean", highlight = "cemean")
  g3 <- plot(ps, type = "indiv")
  g4 <- plot(ps, type = "cusum")

  vdiffr::expect_doppelganger("position-scores", g1)
  vdiffr::expect_doppelganger("position-scores-cemean", g2)
  vdiffr::expect_doppelganger("position-scores-indiv", g3)
  vdiffr::expect_doppelganger("position-scores-cusum", g4)
})
