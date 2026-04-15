test_that("rhs_is_one() works", {
  expect_equal(rhs_is_one(~ 1), TRUE)
  expect_equal(rhs_is_one(X ~ 1), TRUE)
  expect_equal(rhs_is_one(1 ~ X + 1), FALSE)
  expect_equal(rhs_is_one( ~ X + 1), FALSE)
  expect_equal(rhs_is_one(1 ~ X*Y), FALSE)
  expect_equal(rhs_is_one(1 ~ X*Y + 1), FALSE)
})

