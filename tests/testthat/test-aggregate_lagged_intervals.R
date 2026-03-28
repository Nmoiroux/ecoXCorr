test_that("aggregate_lagged_intervals works with a single time series", {
  
  data <- meteoMPL2023
  
  d <- data$date[15:19]
  
  res <- aggregate_lagged_intervals(
    data       = data,
    date_col   = "date",
    value_cols = "rain_sum",
    ref_date   = d,
    interval   = 7,
    max_lag    = 2
  )
  
  expect_true(is.data.frame(res))
  expect_true(all(c("date", "lag_start", "lag_end") %in% names(res)))
  
  # expected number of rows: n_dates * number of lag windows
  n_windows <- sum(sapply(1:2, function(l) 2 - l + 1))
  expect_equal(nrow(res), length(d) * n_windows)
})

test_that("function fails when multiple series are present without id_col", {
  
  data <- meteoMPL2023
  
  # Duplicate dataset → simulate 2 time series
  data2 <- rbind(data, data)
  
  d <- unique(data$date)[15:19]
  
  expect_error(
    aggregate_lagged_intervals(
      data       = data2,
      date_col   = "date",
      value_cols = "rain_sum",
      ref_date   = d,
      interval   = 7,
      max_lag    = 2
    ),
    "Duplicated dates detected"
  )
})

test_that("aggregate_lagged_intervals works with multiple time series", {
  
  data <- meteoMPL2023
  
  # Create artificial multiple series
  data$id <- "A"
  data2 <- data
  data2$id <- "B"
  
  multi_data <- rbind(data, data2)
  
  d <- unique(data$date)[15:19]
  
  res <- aggregate_lagged_intervals(
    data       = multi_data,
    date_col   = "date",
    value_cols = "rain_sum",
    ref_date   = d,
    interval   = 7,
    max_lag    = 2,
    id_col     = "id"
  )
  
  expect_true("id" %in% names(res))
  
  # Check both series are present
  expect_equal(sort(unique(res$id)), c("A", "B"))
  
  # Check row count doubles
  n_windows <- sum(sapply(1:2, function(l) 2 - l + 1))
  expect_equal(nrow(res), 2 * length(d) * n_windows)
})

test_that("results are identical across replicated series", {
  
  data <- meteoMPL2023
  
  # Duplicate identical series
  data$id <- "A"
  data2 <- data
  data2$id <- "B"
  
  multi_data <- rbind(data, data2)
  
  d <- unique(data$date)[15:19]
  
  res <- aggregate_lagged_intervals(
    data       = multi_data,
    date_col   = "date",
    value_cols = "rain_sum",
    ref_date   = d,
    interval   = 7,
    max_lag    = 2,
    id_col     = "id"
  )
  
  res_A <- res[res$id == "A", ]
  res_B <- res[res$id == "B", ]
  
  # Remove id column for comparison
  res_A$id <- NULL
  res_B$id <- NULL
  
  # Remove rownames
  rownames(res_A) <- NULL
  rownames(res_B) <- NULL
  
  expect_equal(res_A, res_B)
})

test_that("aggregated values are correct for a simple known case", {
  
  data <- meteoMPL2023
  
  # pick a small deterministic window
  d <- data$date[15]
  
  res <- aggregate_lagged_intervals(
    data       = data,
    date_col   = "date",
    value_cols = "rain_sum",
    ref_date   = d,
    interval   = 1,
    max_lag    = 1,
    na.rm = TRUE
  )
  
  # expected: previous day
  expected <- data$rain_sum[data$date == (d - 1)]
  
  expect_equal(res$rain_sum_mean[1], expected)
})

test_that("interval length equals interval parameter (7d, 2 lags)", {
  
  data <- meteoMPL2023
  d <- data$date[20]
  
  res <- aggregate_lagged_intervals(
    data       = data,
    date_col   = "date",
    value_cols = "rain_sum",
    ref_date   = d,
    interval   = 7,
    max_lag    = 2
  )
  # we add 1 as end_date (=ref_date when shift = 0) is included
  expect_equal(as.numeric(res$end - res$start + 1), c(7,14,7))
})

test_that("interval length equals interval parameter (1d, 30 lags)", {
  
  data <- meteoMPL2023
  d <- data$date[60]
  
  res <- aggregate_lagged_intervals(
    data       = data,
    date_col   = "date",
    value_cols = "rain_sum",
    ref_date   = d,
    interval   = 1,
    max_lag    = 30
  )
  # we add 1 as end_date (=ref_date when shift = 0) is included
  expect_equal(as.numeric(res$end - res$start + 1)[1:10], 1:10)
})

test_that("aggregate_lagged_intervals returns same values that simple code", {
  
  sampling_dates <- unique(albopictusMPL2023$date)[2:3]
  
  met_agg <- aggregate_lagged_intervals(
    data       = meteoMPL2023,
    date_col   = "date",
    value_cols = c("rain_sum", "temp_mean"),
    ref_date   = sampling_dates,
    interval   = 7,               
    max_lag    = 2,
    shift = 1
  )
  
  shift <- 1
  end <- which(meteoMPL2023$date==sampling_dates[1]) - shift
  start <- which(meteoMPL2023$date==sampling_dates[1]-7) - shift +1
  filter <- meteoMPL2023[start:end,]

  A <- data.frame(
    start = min(filter$date),
    end   = max(filter$date),
    rain_sum_mean = mean(filter$rain_sum),
    rain_sum_min  = min(filter$rain_sum),
    rain_sum_max  = max(filter$rain_sum),
    rain_sum_sum  = sum(filter$rain_sum)
  )
  
  B <- met_agg[1,c(2,3,6:9)]
  
  testthat::expect_equal(A,B)
  
})