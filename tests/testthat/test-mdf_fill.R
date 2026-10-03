# mdf_fill: flip enclosed background to foreground; leave edge-connected
# background and outer patch boundaries untouched.

.fill_mosaik <- function(mat) {
  mosaik(extent = c(0, ncol(mat), 0, nrow(mat)), res = 1, vals = mat)
}

test_that("an enclosed hole is filled", {
  mat <- matrix(c(1,1,1,1,1,
                  1,0,0,0,1,
                  1,0,0,0,1,
                  1,0,0,0,1,
                  1,1,1,1,1), nrow = 5, byrow = TRUE)
  r <- mdf_fill(.fill_mosaik(mat))
  expect_true(all(msk_pull(r, "values") == 1))
})

test_that("edge-connected background is not filled", {
  mat <- matrix(c(1,1,0,1,1,
                  1,1,0,1,1,
                  1,1,1,1,1,
                  1,1,1,1,1,
                  1,1,1,1,1), nrow = 5, byrow = TRUE)
  out <- matrix(msk_pull(mdf_fill(.fill_mosaik(mat)), "values"), 5, 5, byrow = TRUE)
  expect_equal(out[1, 3], 0)   # the notch reaches the top edge -> stays background
  expect_equal(out[2, 3], 0)
})

test_that("the outer boundary of a patch is left intact (not grown)", {
  # a solid 3x3 block of 1s in a 5x5 field of 0s: no enclosed holes, nothing fills
  mat <- matrix(0, 5, 5)
  mat[2:4, 2:4] <- 1
  r <- mdf_fill(.fill_mosaik(mat))
  expect_equal(msk_pull(r, "values"), as.numeric(t(mat)))
})

test_that("holes of any size fill (not just kernel-sized ones)", {
  # a large enclosed hole a morphological close would miss
  mat <- matrix(1, 7, 7)
  mat[2:6, 2:6] <- 0
  mat[1, ] <- 1; mat[7, ] <- 1; mat[, 1] <- 1; mat[, 7] <- 1
  r <- mdf_fill(.fill_mosaik(mat))
  expect_true(all(msk_pull(r, "values") == 1))
})

test_that("a layer that is not binary is refused", {
  mat <- matrix(c(5,5,5,
                  5,0,5,
                  5,5,5), nrow = 3, byrow = TRUE)
  expect_error(mdf_fill(.fill_mosaik(mat)), "not binary")
})

test_that("NA cells count as background and fill when enclosed", {
  mat <- matrix(1, 5, 5)
  mat[3, 3] <- NA
  r <- mdf_fill(.fill_mosaik(mat))
  expect_false(anyNA(msk_pull(r, "values")))
})

test_that("connectivity controls diagonally-pinched holes", {
  # a 2x2 interior hole that touches an edge background cell (top-right corner)
  # ONLY diagonally. Under rook (4) connectivity the diagonal does not connect,
  # so the hole is sealed and fills; under queen (8) it leaks out to the edge
  # and is left as background. The corner cell itself sits on the edge and stays
  # background either way.
  mat <- matrix(c(1,1,1,0,
                  1,0,0,1,
                  1,0,0,1,
                  1,1,1,1), nrow = 4, byrow = TRUE)
  r4 <- matrix(msk_pull(mdf_fill(.fill_mosaik(mat), connectivity = 4L), "values"),
               4, 4, byrow = TRUE)
  r8 <- matrix(msk_pull(mdf_fill(.fill_mosaik(mat), connectivity = 8L), "values"),
               4, 4, byrow = TRUE)

  expect_true(all(r4[2:3, 2:3] == 1))    # trapped under rook -> filled
  expect_true(all(r8[2:3, 2:3] == 0))    # escapes under queen -> not filled
  expect_equal(r4[1, 4], 0)              # edge corner stays background either way
  expect_equal(r8[1, 4], 0)
})

test_that("add writes a new layer and preserves the original", {
  mat <- matrix(c(1,1,1,1,1,
                  1,0,0,0,1,
                  1,0,0,0,1,
                  1,0,0,0,1,
                  1,1,1,1,1), nrow = 5, byrow = TRUE)
  m <- .fill_mosaik(mat)
  r <- mdf_fill(m, add = "filled")
  expect_true("filled" %in% names(r@layers))
  expect_equal(msk_pull(r, "values"), as.numeric(t(mat)))  # original untouched
  expect_true(all(msk_pull(r, "filled") == 1))
})
