test_that("syn_pattern returns a mosaik with correct dimensions", {
  for (type in c("checkerboard", "stripes", "rings", "waves", "grid",
                 "hexagonal")) {
    r <- syn_pattern(mosaik(extent = c(0, 50, 0, 50), res = 1),
                      type = type, frequency = 10)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(50L, 50L))
    expect_equal(length(msk_pull(r)), 2500L)
  }
})

test_that("syn_pattern checkerboard has correct structure", {
  r <- syn_pattern(mosaik(extent = c(0, 20, 0, 20), res = 1),
                    type = "checkerboard", frequency = 10)
  vals <- msk_pull(r)
  expect_true(all(vals %in% c(0, 1)))
  expect_equal(sum(vals == 0), sum(vals == 1))
})

test_that("syn_pattern stripes alternates bands", {
  r <- syn_pattern(mosaik(extent = c(0, 20, 0, 10), res = 1),
                    type = "stripes", frequency = 10, angle = 0)
  vals <- msk_pull(r)
  expect_true(all(vals %in% c(0, 1)))
  mat <- matrix(vals, nrow = 10, ncol = 20, byrow = TRUE)
  expect_equal(length(unique(mat[, 1])), 1)
  expect_equal(length(unique(mat[, 11])), 1)
})

test_that("syn_pattern rings produces continuous values in [0, 1]", {
  r <- syn_pattern(mosaik(extent = c(0, 50, 0, 50), res = 1),
                    type = "rings", frequency = 3)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  expect_equal(mat[25, 25], 1, tolerance = 0.1)
})

test_that("syn_pattern waves produces continuous values in [0, 1]", {
  r <- syn_pattern(mosaik(extent = c(0, 50, 0, 50), res = 1),
                    type = "waves", frequency = 4)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_pattern grid marks edges correctly", {
  r <- syn_pattern(mosaik(extent = c(0, 20, 0, 20), res = 1),
                    type = "grid", frequency = 10)
  vals <- msk_pull(r)
  expect_true(all(vals %in% c(0, 1)))
  expect_equal(vals[1], 1)
  mat <- matrix(vals, nrow = 20, ncol = 20, byrow = TRUE)
  expect_equal(mat[6, 6], 0)
})

test_that("syn_pattern hexagonal produces valid tile IDs", {
  r <- syn_pattern(mosaik(extent = c(0, 50, 0, 50), res = 1),
                    type = "hexagonal", frequency = 8)
  vals <- msk_pull(r)
  expect_true(all(vals %in% c(0, 1, 2)))
})

test_that("syn_pattern name sets layer name", {
  r <- syn_pattern(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    type = "checkerboard", name = "pattern")
  expect_equal(names(r@layers), "pattern")
})

test_that("syn_pattern records provenance", {
  r <- syn_pattern(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    type = "rings", frequency = 2)
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "syn_pattern")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "rings")
})
