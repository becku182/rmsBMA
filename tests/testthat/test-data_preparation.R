test_that("standardization only (no ids) produces mean 0, sd 1 for numeric vars", {
  set.seed(1)
  df <- data.frame(
    a = rnorm(50, mean = 10, sd = 2),
    b = rnorm(50, mean = -3, sd = 5),
    c = letters[1:50]  # non-numeric
  )

  X <- data_preparation(df, standardize = TRUE)

  expect_true(is.data.frame(X))
  expect_setequal(names(X), c("a", "b", "c"))  # non-numeric untouched/kept

  expect_equal(mean(X$a, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(sd(X$a, na.rm = TRUE),   1, tolerance = 1e-12)

  expect_equal(mean(X$b, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(sd(X$b, na.rm = TRUE),   1, tolerance = 1e-12)
})

test_that("no ids + standardize = FALSE errors", {
  df <- data.frame(a = 1:5)
  expect_error(
    data_preparation(df, standardize = FALSE),
    "set `standardize = TRUE`"
  )
})

test_that("errors when no numeric variables exist", {
  df <- data.frame(x = letters[1:5], y = LETTERS[1:5])
  expect_error(data_preparation(df, standardize = TRUE), "No numeric variables")
})

test_that("with ids present: ids are dropped and only non-id numeric vars transformed", {
  set.seed(2)
  df <- data.frame(
    id = rep(1:5, each = 4),
    time = rep(1:4, times = 5),
    y = rnorm(20),
    z = rnorm(20),
    w = letters[1:20]
  )

  X <- data_preparation(df, id = "id", time = "time",
                        fixed_effects = FALSE, standardize = TRUE)

  expect_false("id" %in% names(X))
  expect_false("time" %in% names(X))
  expect_true(all(c("y", "z", "w") %in% names(X)))

  # y, z standardized
  expect_equal(mean(X$y, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(sd(X$y, na.rm = TRUE),   1, tolerance = 1e-12)

  expect_equal(mean(X$z, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(sd(X$z, na.rm = TRUE),   1, tolerance = 1e-12)

  # non-numeric unchanged
  expect_identical(X$w, df$w)
})

test_that("section FE demeaning yields zero group means by id (for transformed vars)", {
  set.seed(3)
  df <- data.frame(
    id = rep(1:4, each = 5),
    time = rep(1:5, times = 4),
    x = rnorm(20, mean = 2),
    y = rnorm(20, mean = -1)
  )

  X <- data_preparation(df, id = "id", time = "time",
                        fixed_effects = TRUE, effect = "section",
                        standardize = FALSE)

  # group means by id should be ~0 after section FE
  mx <- tapply(X$x, df$id, mean, na.rm = TRUE)
  my <- tapply(X$y, df$id, mean, na.rm = TRUE)

  expect_true(all(abs(mx) < 1e-12))
  expect_true(all(abs(my) < 1e-12))
})

test_that("time FE demeaning yields zero group means by time (for transformed vars)", {
  set.seed(4)
  df <- data.frame(
    id = rep(1:5, each = 4),
    time = rep(1:4, times = 5),
    x = rnorm(20),
    y = rnorm(20)
  )

  X <- data_preparation(df, id = "id", time = "time",
                        fixed_effects = TRUE, effect = "time",
                        standardize = FALSE)

  mx <- tapply(X$x, df$time, mean, na.rm = TRUE)
  my <- tapply(X$y, df$time, mean, na.rm = TRUE)

  expect_true(all(abs(mx) < 1e-12))
  expect_true(all(abs(my) < 1e-12))
})

test_that("two-way FE demeaning yields zero means by id and by time (balanced panel)", {
  set.seed(5)
  df <- data.frame(
    id = rep(1:6, each = 4),
    time = rep(1:4, times = 6),
    x = rnorm(24, mean = 10),
    y = rnorm(24, mean = -2)
  )

  X <- data_preparation(df, id = "id", time = "time",
                        fixed_effects = TRUE, effect = "twoway",
                        standardize = FALSE)

  mx_id <- tapply(X$x, df$id, mean, na.rm = TRUE)
  mx_t  <- tapply(X$x, df$time, mean, na.rm = TRUE)

  my_id <- tapply(X$y, df$id, mean, na.rm = TRUE)
  my_t  <- tapply(X$y, df$time, mean, na.rm = TRUE)

  expect_true(all(abs(mx_id) < 1e-12))
  expect_true(all(abs(mx_t)  < 1e-12))
  expect_true(all(abs(my_id) < 1e-12))
  expect_true(all(abs(my_t)  < 1e-12))
})

test_that("FE + standardization gives (approx) mean 0 and sd 1", {
  set.seed(6)
  df <- data.frame(
    id = rep(1:8, each = 3),
    time = rep(1:3, times = 8),
    x = rnorm(24, mean = 7, sd = 4),
    y = rnorm(24, mean = -5, sd = 2)
  )

  X <- data_preparation(df, id = "id", time = "time",
                        fixed_effects = TRUE, effect = "twoway",
                        standardize = TRUE)

  expect_equal(mean(X$x, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(sd(X$x, na.rm = TRUE),   1, tolerance = 1e-12)

  expect_equal(mean(X$y, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(sd(X$y, na.rm = TRUE),   1, tolerance = 1e-12)
})

test_that("standardization sets constant / zero-variance vars to NA_real_", {
  df <- data.frame(
    id = rep(1:3, each = 4),
    time = rep(1:4, times = 3),
    x = rep(5, 12),              # constant
    y = c(1:11, NA_real_)        # non-constant with NA
  )

  X <- data_preparation(df, id = "id", time = "time",
                        fixed_effects = FALSE, standardize = TRUE)

  expect_true(all(is.na(X$x)))
  # y should be standardized (not all NA)
  expect_true(any(!is.na(X$y)))
  expect_equal(mean(X$y, na.rm = TRUE), 0, tolerance = 1e-12)
})

test_that("errors if ids provided but not found in data", {
  df <- data.frame(id = 1:5, time = 1:5, x = rnorm(5))
  expect_error(data_preparation(df, id = "ID", time = "time", standardize = TRUE))
  expect_error(data_preparation(df, id = "id", time = "TIME", standardize = TRUE))
})

test_that("errors if only id or only time is provided and standardize = FALSE", {
  df <- data.frame(id = 1:5, time = 1:5, x = rnorm(5))

  expect_error(
    data_preparation(df, id = "id", standardize = FALSE),
    "set `standardize = TRUE`"
  )

  expect_error(
    data_preparation(df, time = "time", standardize = FALSE),
    "set `standardize = TRUE`"
  )
})

test_that("with ids: errors when no numeric vars besides id/time", {
  df <- data.frame(
    id = 1:5,
    time = 1:5,
    group = letters[1:5]
  )

  expect_error(
    data_preparation(df, id = "id", time = "time", standardize = TRUE),
    "No numeric variables to transform"
  )
})
