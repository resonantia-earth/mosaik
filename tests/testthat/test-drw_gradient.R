test_that("drw_gradient returns a mosaik with correct dimensions", {
  for (type in c("planar", "point", "line", "circle", "square")) {
    r <- drw_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                       type = type)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(50L, 50L))
    expect_equal(length(msk_pull(r)), 2500L)
  }
})

test_that("drw_gradient planar values are in [0, 1]", {
  r <- drw_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "planar")
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_true(mat[25, 1] < 0.1)
  expect_true(mat[25, 50] > 0.9)
})

test_that("drw_gradient planar angle rotates correctly", {
  r <- drw_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "planar", angle = 90)
  vals <- msk_pull(r)
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  # 90 degrees increases from bottom to top; row 1 of the matrix is the top
  expect_true(mean(mat[1, ]) > 0.9)
  expect_true(mean(mat[50, ]) < 0.1)
})

test_that("drw_gradient position counts from the lower-left corner", {
  r <- drw_gradient(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    type = "point", position = c(0.25, 0.25))
  v <- msk_pull(r)
  # the minimum lies in the lower left
  expect_equal(which.min(v), .cell(r, x = 2.5, y = 2.5))
})

test_that("drw_gradient point peaks at centre", {
  r <- drw_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "point", invert = TRUE)
  vals <- msk_pull(r)
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_true(mat[25, 25] > 0.9)
  expect_true(mat[1, 1] < 0.2)
})

test_that("drw_gradient line creates symmetric distance field", {
  r <- drw_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "line", angle = 0)
  vals <- msk_pull(r)
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_equal(mat[1, 25], mat[50, 25], tolerance = 0.05)
})

test_that("drw_gradient circle creates radial distance from ring", {
  r <- drw_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "circle", size = 0.3)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("drw_gradient invert works", {
  r1 <- drw_gradient(mosaik(extent = c(0, 30, 0, 30), res = 1),
                      type = "planar")
  r2 <- drw_gradient(mosaik(extent = c(0, 30, 0, 30), res = 1),
                      type = "planar", invert = TRUE)
  v1 <- msk_pull(r1)
  v2 <- msk_pull(r2)
  expect_equal(v1 + v2, rep(1, length(v1)), tolerance = 1e-10)
})

test_that("drw_gradient takes its origin from a layer of the object", {
  origin_vals <- rep(0, 400)
  origin_vals[190] <- 1
  m <- mosaik(extent = c(0, 20, 0, 20), res = 1,
              vals = list(v = origin_vals))
  r <- drw_gradient(m, origin = "v", name = "g")
  vals <- msk_pull(r, "g")
  expect_true(all(vals >= 0 & vals <= 1))
  expect_equal(vals[190], 0)
  expect_error(drw_gradient(m, origin = "nothere"), "names no layer")
  expect_error(drw_gradient(msk_add(m, w = rep(2, 400)), origin = "w"),
               "binary")
})

test_that("drw_gradient turns the square by angle", {
  m <- mosaik(extent = c(0, 41, 0, 41), res = 1)
  s0 <- msk_pull(drw_gradient(m, type = "square", size = 0.5))
  s45 <- msk_pull(drw_gradient(m, type = "square", size = 0.5, angle = 45))
  # the corner of the upright square is inside it, the turned one's is not
  corner <- .cell(m, x = 30, y = 30)
  expect_equal(s0[corner], 0)
  expect_gt(s45[corner], 0)
})

test_that("drw_gradient name sets layer name", {
  r <- drw_gradient(mosaik(extent = c(0, 10, 0, 10), res = 1),
                     type = "planar", name = "gradient")
  expect_equal(names(r@layers), "gradient")
})

test_that("drw_gradient records provenance", {
  r <- drw_gradient(mosaik(extent = c(0, 10, 0, 10), res = 1),
                     type = "point")
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "drw_gradient")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "point")
})
