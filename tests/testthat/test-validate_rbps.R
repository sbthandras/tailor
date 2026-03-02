data(rbps)

test_that("validate_rbps", {
  expect_true(validate_rbps(rbps))

  rbps2 <- rbps
  names(rbps2)[which(names(rbps2) == "Core_ORF")] <- "noname"
  expect_false(validate_rbps(rbps2))
  expect_warning(validate_rbps(rbps2, verbose = TRUE),
                 "Variable 'Core_ORF' not found.")

  rbps2 <- rbps
  rbps2$Core_ORF[2] = rbps2$Core_ORF[1]
  expect_false(validate_rbps(rbps2))
  expect_warning(validate_rbps(rbps2, verbose = TRUE),
                 "Variable 'Core_ORF' not unique.")

  rbps2 <- rbps
  names(rbps2)[which(names(rbps2) == "translation")] <- "noname"
  expect_false(validate_rbps(rbps2))
  expect_warning(validate_rbps(rbps2, verbose = TRUE),
                 "Variable 'translation' not found.")

  rbps2 <- rbps
  rbps2$translation <- seq_along(rbps2$translation)
  expect_false(validate_rbps(rbps2))
  expect_warning(validate_rbps(rbps2, verbose = TRUE),
                 "Variable 'translation' is not character.")

  rbps2 <- rbps
  rbps2$translation[1:3] <- ""
  expect_false(validate_rbps(rbps2))
  expect_warning(validate_rbps(rbps2, verbose = TRUE),
                 paste0("Variable 'translation' contains an empty string in ",
                        "the following rows: 1, 2, 3"))

  rbps2 <- rbps
  rbps2$translation[1:3] <- paste0(rbps2$translation[1:3], "*")
  expect_false(validate_rbps(rbps2))

  expr <- rlang::expr(validate_rbps(rbps2, verbose = TRUE))
  expect_warning(eval(expr))

  msg <- capture_warning(eval(expr))
  expect_equal(
    msg$message,
    paste0(
      "Variable 'translation' contains unknown ",
      "character(s) in the following rows: 1, 2, 3"
    )
  )

})
