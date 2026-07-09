test_that("Breslow baseline evaluation handles empty hazards", {
  bh <- data.frame(hazard = numeric(), time = numeric())

  H0 <- .basehaz_at_horizons(bh, c(1, 3, 5))

  expect_equal(H0, c(0, 0, 0))
})

test_that("Breslow baseline evaluation handles a single hazard point", {
  bh <- data.frame(hazard = 0.25, time = 2)

  H0 <- .basehaz_at_horizons(bh, c(1, 2, 5))

  expect_equal(H0, c(0, 0.25, 0.25))
})

test_that("weighted KM evaluation handles no-event samples", {
  S_t <- .km_at(time = c(1, 2, 3), status = c(0, 0, 0), t = 2)

  expect_equal(S_t, 1)
})
