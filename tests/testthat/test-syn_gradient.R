test_that("syn_gradient returns a mosaik with correct dimensions", {
  for (type in c("planar", "point", "line", "circle", "square", "triangle",
                 "hexagon")) {
    r <- syn_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                       type = type)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(50L, 50L))
    expect_equal(length(msk_pull(r)), 2500L)
  }
})

test_that("syn_gradient planar values are in [0, 1]", {
  r <- syn_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "planar")
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_true(mat[25, 1] < 0.1)
  expect_true(mat[25, 50] > 0.9)
})

test_that("syn_gradient planar angle rotates correctly", {
  r <- syn_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "planar", angle = 90)
  vals <- msk_pull(r)
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_true(mean(mat[1, ]) < 0.1)
  expect_true(mean(mat[50, ]) > 0.9)
})

test_that("syn_gradient point peaks at centre", {
  r <- syn_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "point", invert = TRUE)
  vals <- msk_pull(r)
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_true(mat[25, 25] > 0.9)
  expect_true(mat[1, 1] < 0.2)
})

test_that("syn_gradient line creates symmetric distance field", {
  r <- syn_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "line", angle = 0)
  vals <- msk_pull(r)
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_equal(mat[1, 25], mat[50, 25], tolerance = 0.05)
})

test_that("syn_gradient circle creates radial distance from ring", {
  r <- syn_gradient(mosaik(extent = c(0, 50, 0, 50), res = 1),
                     type = "circle", size = 0.3)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_gradient invert works", {
  r1 <- syn_gradient(mosaik(extent = c(0, 30, 0, 30), res = 1),
                      type = "planar")
  r2 <- syn_gradient(mosaik(extent = c(0, 30, 0, 30), res = 1),
                      type = "planar", invert = TRUE)
  v1 <- msk_pull(r1)
  v2 <- msk_pull(r2)
  expect_equal(v1 + v2, rep(1, length(v1)), tolerance = 1e-10)
})

test_that("syn_gradient origin mosaik works", {
  origin_vals <- rep(0, 400)
  origin_vals[190] <- 1
  origin <- mosaik(extent = c(0, 20, 0, 20), res = 1,
                   vals = list(v = origin_vals))
  r <- syn_gradient(mosaik(extent = c(0, 20, 0, 20), res = 1),
                     origin = origin)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
  expect_equal(vals[190], 0)
})

test_that("syn_gradient name sets layer name", {
  r <- syn_gradient(mosaik(extent = c(0, 10, 0, 10), res = 1),
                     type = "planar", name = "gradient")
  expect_equal(names(r@layers), "gradient")
})

test_that("syn_gradient records provenance", {
  r <- syn_gradient(mosaik(extent = c(0, 10, 0, 10), res = 1),
                     type = "point")
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "syn_gradient")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "point")
})
