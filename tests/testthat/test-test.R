test_that("aggregate_lagged_intervals returns expected number of lags", {
  i <- sample(c(1:10), 1)
  m <- sample(c(1:10), 1)
  x <-sample(c(1:5), 1)

  res <- aggregate_lagged_intervals(
    data = meteoMPL2023,
    date_col = "date",
    value_cols = "rain_sum",
    ref_date = unique(albopictusMPL2023$date),
    interval = i,
    max_lag = m,
    shift = x,
    funs = list(sum = sum)
  )

  expect_true(all(c("lag_start", "lag_end") %in% names(res)))
  expect_equal(as.integer(max(res$end-res$start))+1, i*m)
  expect_equal(length(unique(albopictusMPL2023$date)) * m * (m+1)/2, nrow(res))
})

