grid2 <- function(){
  a <- matrix(c(0,1,1, 1,0,1, 1,1,0), 3, 3, byrow = TRUE)
  b <- matrix(c(0,0,1, 1,0,0, 0,1,1), 3, 3, byrow = TRUE)
  mosaik(extent = c(0, 3, 0, 3), res = 1,
         vals = list(a = as.numeric(t(a)), b = as.numeric(t(b))))
}

test_that("value = TRUE keeps the layer's own values and NAs the rest", {
  m <- grid2()
  o <- mdf_filter(m, a == 1 & b == 0, value = TRUE, layer = "a")
  expect_equal(msk_pull(o, "a"), c(NA,1,NA, NA,NA,1, 1,NA,NA))
})

test_that("the default writes the predicate as a 1/0 mask", {
  m <- grid2()
  # difference, union and intersection are all the same operator
  d <- mdf_filter(m, a == 1 & b == 0, add = "d")
  u <- mdf_filter(m, a == 1 | b == 1, add = "u")
  i <- mdf_filter(m, a == 1 & b == 1, add = "i")

  expect_equal(msk_pull(d, "d"), c(0,1,0, 0,0,1, 1,0,0))
  expect_equal(msk_pull(u, "u"), c(0,1,1, 1,0,1, 1,1,1))
  expect_equal(msk_pull(i, "i"), c(0,0,1, 1,0,0, 0,1,0))
})

test_that("value takes only TRUE or FALSE", {
  m <- grid2()
  expect_error(mdf_filter(m, a == 1, value = 7, add = "seven"))
})

test_that("a mask does not read the destination layer's values", {
  m <- grid2()
  # filter on b but write into a new layer: a must be untouched, and the
  # result must be b's mask rather than a's values
  o <- mdf_filter(m, b == 1, layer = "a", add = "fromB")
  expect_equal(msk_pull(o, "fromB"), c(0,0,1, 1,0,0, 0,1,1))
  expect_equal(msk_pull(o, "a"), msk_pull(m, "a"))
})

test_that("add = NULL overwrites the layer in place", {
  m <- grid2()
  o <- mdf_filter(m, a == 1, layer = "a")
  expect_equal(msk_pull(o, "a"), c(0,1,1, 1,0,1, 1,1,0))
  expect_equal(names(o@layers), names(m@layers))
})

test_that("an always-false predicate gives an empty mask", {
  m <- grid2()
  o <- mdf_filter(m, a < 0, add = "none")
  expect_equal(sum(msk_pull(o, "none")), 0)
})

test_that("a non-logical predicate is rejected", {
  m <- grid2()
  expect_error(mdf_filter(m, a + b), "logical")
})

test_that("a predicate is recorded unevaluated and replays correctly", {
  # .record_step evaluates arguments so a step stores objects rather than
  # calls, but a predicate names layers of the mosaik the recipe is applied to
  # LATER, so it must survive recording unevaluated.
  r <- mdf_filter(obj = NULL, a == 1 & b == 0, add = "d")
  expect_equal(deparse(r@provenance[[1]][[1]]$wasGeneratedBy$withArguments$expr), "a == 1 & b == 0")

  m <- grid2()
  expect_equal(msk_pull(mdf(m, r), "d"), c(0,1,0, 0,0,1, 1,0,0))
})

test_that("a recipe mixes quoted predicates with evaluated arguments", {
  r <- mdf_dilate(obj = NULL, struct = msk_struct("square", width = 3),
                  layer = "a", add = "dil") |>
       mdf_filter(dil == 1 & b == 0, add = "res")

  # struct is stored as the built object, expr as an unevaluated expression
  expect_s4_class(r@provenance[[1]][[1]]$wasGeneratedBy$withArguments$struct, "struct")
  expect_true(is.call(r@provenance[[2]][[1]]$wasGeneratedBy$withArguments$expr))

  m <- grid2()
  expect_equal(msk_pull(mdf(m, r), "res"), c(1,1,0, 0,1,1, 1,0,0))
})
