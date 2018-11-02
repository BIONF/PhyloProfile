context("test tmp")

test_that("just for testing testthat", {
  a <- 2
  expect_that( a, is_more_than(0))
})
