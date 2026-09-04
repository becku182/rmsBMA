test_that("group_dilution returns 1 for all models when there are no groups", {
  Reg_ID <- rbind(
    c(0, 0, 0),
    c(1, 0, 1),
    c(1, 1, 1)
  )
  Nar_vec <- c(0, 0, 0)

  out <- group_dilution(Reg_ID, Nar_vec, p = 0.7)

  expect_type(out, "double")
  expect_length(out, nrow(Reg_ID))
  expect_equal(out, rep(1, nrow(Reg_ID)))
})

test_that("group_dilution works with scalar p (same for all groups)", {
  # K=5, groups: 1 for cols 1:3, 2 for cols 4:5
  Nar_vec <- c(1, 1, 1, 2, 2)

  Reg_ID <- rbind(
    c(0,0,0,0,0),  # null: D=0 -> 1
    c(1,0,0,0,0),  # one from g1: D=0 -> 1
    c(1,1,0,0,0),  # two from g1: D=1 -> p
    c(1,1,1,0,0),  # three from g1: D=2 -> p^2
    c(1,1,0,1,1)   # g1 count=2 -> +1; g2 count=2 -> +1 => D=2 -> p^2
  )

  p <- 0.7
  out <- group_dilution(Reg_ID, Nar_vec, p = p)

  expected <- c(
    1,
    1,
    p^1,
    p^2,
    p^2
  )

  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("group_dilution works with group-specific p vector in group order", {
  Nar_vec <- c(1, 1, 1, 2, 2)
  Reg_ID <- rbind(
    c(1,1,0,0,0),  # g1 count=2 -> expo1=1 => p1
    c(1,1,1,0,0),  # g1 count=3 -> expo1=2 => p1^2
    c(0,0,0,1,1),  # g2 count=2 -> expo2=1 => p2
    c(1,1,0,1,1)   # expo1=1 and expo2=1 => p1*p2
  )

  p1 <- 0.7
  p2 <- 0.5
  out <- group_dilution(Reg_ID, Nar_vec, p = c(p1, p2))

  expected <- c(
    p1,
    p1^2,
    p2,
    p1 * p2
  )

  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("group_dilution matches named p to group IDs", {
  Nar_vec <- c(2, 2, 1, 1)  # groups are {1,2} after sorting
  Reg_ID <- rbind(
    c(1,1,0,0),  # group 2 count=2 -> expo2=1
    c(0,0,1,1),  # group 1 count=2 -> expo1=1
    c(1,1,1,1)   # group2 expo1 + group1 expo1 => p2 * p1
  )

  p_named <- c("1" = 0.6, "2" = 0.2)
  out <- group_dilution(Reg_ID, Nar_vec, p = p_named)

  expected <- c(
    0.2,        # group 2 penalty
    0.6,        # group 1 penalty
    0.2 * 0.6
  )

  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("group_dilution handles p = 0 (penalize repeated group members to 0)", {
  Nar_vec <- c(1, 1, 2)

  Reg_ID <- rbind(
    c(1,0,1),  # group1 count=1 => expo=0 => no penalty => 1
    c(1,1,0),  # group1 count=2 => expo=1 => p^1 => 0
    c(1,1,1)   # group1 expo=1 => 0 regardless of group2
  )

  out <- group_dilution(Reg_ID, Nar_vec, p = 0)

  expect_equal(out, c(1, 0, 0))
})

test_that("group_dilution validates inputs", {
  Reg_ID <- rbind(c(0,1,0), c(1,1,0))
  Nar_vec <- c(1, 1, 2)

  # Nar_vec length mismatch
  expect_error(group_dilution(Reg_ID, Nar_vec[-1], p = 0.5), "Nar_vec must have length")

  # p must be numeric and not NA
  expect_error(group_dilution(Reg_ID, Nar_vec, p = NA_real_), "p must be numeric")
  expect_error(group_dilution(Reg_ID, Nar_vec, p = "0.5"), "p must be numeric")

  # p in [0,1]
  expect_error(group_dilution(Reg_ID, Nar_vec, p = -0.1), "\\[0, 1\\]")
  expect_error(group_dilution(Reg_ID, Nar_vec, p =  1.1), "\\[0, 1\\]")

  # p length mismatch: groups are {1,2} => G=2
  expect_error(group_dilution(Reg_ID, Nar_vec, p = c(0.5, 0.4, 0.3)), "length 1 or length equal")

  # named p missing group ids
  expect_error(group_dilution(Reg_ID, Nar_vec, p = c("1" = 0.5, "3" = 0.2)),
               "not all group IDs")
})
