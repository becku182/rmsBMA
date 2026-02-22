test_that("subset_design selects columns where model_row != 0", {
  x <- matrix(1:20, nrow = 5, ncol = 4)
  colnames(x) <- c("a", "b", "c", "d")

  model_row <- c(1, 0, 2, 0)  # nonzero means included

  Xsub <- subset_design(x, model_row)

  expect_true(is.matrix(Xsub))
  expect_equal(Xsub, x[, c(1, 3), drop = FALSE])
  expect_equal(colnames(Xsub), c("a", "c"))
})

test_that("subset_design preserves column order and does not drop to vector", {
  set.seed(1)
  x <- matrix(rnorm(30), nrow = 10, ncol = 3)

  model_row <- c(0, 1, 0)
  Xsub <- subset_design(x, model_row)

  expect_true(is.matrix(Xsub))
  expect_equal(dim(Xsub), c(nrow(x), 1))
  expect_equal(Xsub, x[, 2, drop = FALSE])
})

test_that("subset_design returns a 0-column matrix if nothing selected", {
  set.seed(2)
  x <- matrix(rnorm(24), nrow = 6, ncol = 4)

  model_row <- c(0, 0, 0, 0)
  Xsub <- subset_design(x, model_row)

  expect_true(is.matrix(Xsub))
  expect_equal(nrow(Xsub), nrow(x))
  expect_equal(ncol(Xsub), 0)
})

test_that("subset_design works with logical inclusion vectors", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- c("a", "b", "c", "d")

  model_row <- c(TRUE, FALSE, TRUE, FALSE)
  Xsub <- subset_design(x, model_row)

  expect_equal(Xsub, x[, c(1, 3), drop = FALSE])
  expect_equal(colnames(Xsub), c("a", "c"))
})

test_that("subset_design with shorter model_row ignores columns beyond its length", {
  # Current behavior: no error; it can only select among 1:length(model_row)
  set.seed(3)
  x <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colnames(x) <- c("a", "b", "c", "d")

  model_row <- c(1, 0, 1)  # length 3, but x has 4 columns
  Xsub <- subset_design(x, model_row)

  expect_equal(Xsub, x[, c(1, 3), drop = FALSE])
  expect_equal(colnames(Xsub), c("a", "c"))
})

test_that("subset_design errors when model_row refers to columns that don't exist (longer vector)", {
  # If model_row is longer than ncol(x) and has nonzero entries beyond ncol(x),
  # idx will contain out-of-range column indices -> subsetting error.
  set.seed(4)
  x <- matrix(rnorm(20), nrow = 5, ncol = 4)

  model_row <- c(0, 1, 0, 1, 1)  # includes column 5, which doesn't exist
  expect_error(subset_design(x, model_row))
})
