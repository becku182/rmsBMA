test_that("model_matrix returns a matrix with correct dimensions", {
  K <- 5
  M <- 3

  mm <- model_matrix(K, M)

  expect_true(is.matrix(mm))
  expect_equal(ncol(mm), K)
  expect_equal(nrow(mm), sum(choose(K, 0:M)))
})

test_that("model_matrix contains the null model as the first row", {
  K <- 6
  M <- 2

  mm <- model_matrix(K, M)

  expect_equal(mm[1, ], rep(0L, K))
})

test_that("model_matrix uses only 0/1 entries and models have size <= M", {
  K <- 7
  M <- 3

  mm <- model_matrix(K, M)

  expect_true(all(mm %in% c(0L, 1L)))
  expect_true(all(rowSums(mm) <= M))
  expect_true(all(rowSums(mm) >= 0))
})

test_that("model_matrix has the correct number of models of each size 0..M", {
  K <- 6
  M <- 4

  mm <- model_matrix(K, M)
  tab <- table(rowSums(mm))

  # for k = 0..M, should have choose(K,k) models
  for (k in 0:M) {
    expect_equal(as.integer(tab[as.character(k)]), choose(K, k))
  }
})

test_that("model_matrix ordering is: null model, then all 1-regressor models, then 2-regressor models, ...", {
  K <- 4
  M <- 3

  mm <- model_matrix(K, M)

  # first is null
  expect_equal(mm[1, ], c(0, 0, 0, 0))

  # next K rows should be the singletons in order (1,2,3,4)
  expected_singletons <- rbind(
    c(1,0,0,0),
    c(0,1,0,0),
    c(0,0,1,0),
    c(0,0,0,1)
  )
  expect_equal(mm[2:(K+1), ], expected_singletons)
})

test_that("model_matrix exact structure for small case K=3, M=2", {
  K <- 3
  M <- 2

  mm <- model_matrix(K, M)

  expected <- rbind(
    c(0,0,0),
    c(1,0,0),
    c(0,1,0),
    c(0,0,1),
    c(1,1,0),
    c(1,0,1),
    c(0,1,1)
  )

  expect_equal(mm, expected)
})

test_that("model_matrix errors when M > K", {
  expect_error(model_matrix(3, 4))
})
