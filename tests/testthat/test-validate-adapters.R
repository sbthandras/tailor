data(adapters)

test_that("adapter objects remain unchanged", {
  ada <- adapters |> as.data.frame() |> validate_adapters()
  expect_true(all.equal(adapters, ada))
})

test_that("validation removes duplicates", {
  ada <- adapters |>
    as.data.frame() |>
    dplyr::bind_rows(
      data.frame(
        pattern_id = "MN395291-1",
        subject_id = "ON513429-1",
        start = 1,
        end = 177,
        mean_score = 4.323,
        pident = 0.876
      )
    )
  expect_true(all.equal(adapters, ada |> validate_adapters()))
})

test_that("validation removes reverse duplicates", {
  ada <- adapters |>
    as.data.frame() |>
    dplyr::bind_rows(
      data.frame(
        pattern_id = "ON513429-1",
        subject_id = "MN395291-1",
        start = 1,
        end = 177,
        mean_score = 4.323,
        pident = 0.876
      )
    )
  expect_true(all.equal(adapters, ada |> validate_adapters()))
})

test_that("Throw an error if a pair has multiple values - forward", {
  ada <- adapters |>
    as.data.frame() |>
    dplyr::bind_rows(
      data.frame(
        pattern_id = "MN395291-1",
        subject_id = "ON513429-1",
        start = 1,
        end = 177,
        mean_score = 4.323,
        pident = 0.9 # new value
      )
    )
  expect_error(
    ada |> validate_adapters(),
    "Unresolved duplicate pairs detected: MN395291-1 ON513429-1"
  )
})

test_that("Throw and error if a pair has multiple values - reverse", {
  ada <- adapters |>
    as.data.frame() |>
    dplyr::bind_rows(
      data.frame(
        pattern_id = "ON513429-1",
        subject_id = "MN395291-1",
        start = 1,
        end = 177,
        mean_score = 4.323,
        pident = 0.9 # new value
      )
    )
  expect_error(
    ada |> validate_adapters(),
    "Unresolved duplicate pairs detected: MN395291-1 ON513429-1"
  )
})


