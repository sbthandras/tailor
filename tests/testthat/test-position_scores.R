data(rbps)

test_that("position_scores",{
  ps <- position_scores(rbps$Core_ORF[1], rbps$Core_ORF[2], data = rbps)

  expect_s3_class(ps, c("ps", "list"))
  expect_length(ps, 4)
  expect_type(ps$pattern_id, "character")
  expect_type(ps$subject_id, "character")
  expect_s3_class(ps$position_scores, c("tbl_df", "tbl", "data.frame"))
  expect_equal(nrow(ps$position_scores), 913)
})

test_that("position_scores throws an error when sequence is empty", {
  rbps2 <- dplyr::bind_rows(
    rbps,
    data.frame("Core_ORF" = "neworf", "translation" = "")
  )
  expr <- rlang::expr(position_scores(rbps$Core_ORF[1], "neworf", data = rbps2))
  msg <- capture_error(eval(expr))

  expect_error(eval(expr))
  expect_equal(
    msg$message,
    "Argument 'data' is not valid. See validate_rbps() for more information.")
})

test_that("position_scores no longer fails on trailing '*' in sequence", {
  rbps2 <- dplyr::bind_rows(
    rbps,
    data.frame("Core_ORF" = "neworf", "translation" = "AAQ*ANA*")
  )
  ps <- position_scores(rbps$Core_ORF[1], "neworf", data = rbps2)

  expect_s3_class(ps, c("ps", "list"))
  expect_s3_class(ps$position_scores, c("tbl_df", "tbl", "data.frame"))
  expect_type(ps$position_scores$score, "integer")
})
