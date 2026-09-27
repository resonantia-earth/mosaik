test_that("syn_noise returns a mosaik with correct dimensions", {
  for (type in c("white", "perlin", "simplex")) {
    r <- syn_noise(mosaik(extent = c(0, 50, 0, 50), res = 1),
                    type = type, seed = 42)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(50L, 50L))
    expect_equal(length(msk_pull(r)), 2500L)
  }
})

test_that("syn_noise white produces uniform values in [0, 1]", {
  r <- syn_noise(mosaik(extent = c(0, 100, 0, 100), res = 1),
                  type = "white", seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
  expect_equal(mean(vals), 0.5, tolerance = 0.05)
})

test_that("syn_noise perlin produces values in [0, 1]", {
  r <- syn_noise(mosaik(extent = c(0, 100, 0, 100), res = 1),
                  type = "perlin", frequency = 4, seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_noise perlin is spatially autocorrelated", {
  r <- syn_noise(mosaik(extent = c(0, 50, 0, 50), res = 1),
                  type = "perlin", frequency = 4, seed = 1)
  vals <- msk_pull(r)
  mat <- matrix(vals, nrow = 50, ncol = 50, byrow = TRUE)
  h_diff <- mean(abs(diff(mat[25, ])))
  rand_diff <- mean(abs(diff(sample(vals))))
  expect_lt(h_diff, rand_diff)
})

test_that("syn_noise frequency controls scale", {
  r_low <- syn_noise(mosaik(extent = c(0, 100, 0, 100), res = 1),
                      type = "perlin", frequency = 2, seed = 1)
  r_high <- syn_noise(mosaik(extent = c(0, 100, 0, 100), res = 1),
                       type = "perlin", frequency = 16, seed = 1)
  mat_low <- matrix(msk_pull(r_low), nrow = 100, ncol = 100, byrow = TRUE)
  mat_high <- matrix(msk_pull(r_high), nrow = 100, ncol = 100, byrow = TRUE)
  var_low <- mean(abs(diff(mat_low[50, ])))
  var_high <- mean(abs(diff(mat_high[50, ])))
  expect_lt(var_low, var_high)
})

test_that("syn_noise is reproducible with seed", {
  r1 <- syn_noise(mosaik(extent = c(0, 30, 0, 30), res = 1),
                   type = "perlin", seed = 42)
  r2 <- syn_noise(mosaik(extent = c(0, 30, 0, 30), res = 1),
                   type = "perlin", seed = 42)
  expect_equal(msk_pull(r1), msk_pull(r2))
})

test_that("syn_noise name sets layer name", {
  r <- syn_noise(mosaik(extent = c(0, 10, 0, 10), res = 1),
                  type = "white", name = "noise", seed = 1)
  expect_equal(names(r@layers), "noise")
})

test_that("syn_noise records provenance", {
  r <- syn_noise(mosaik(extent = c(0, 10, 0, 10), res = 1),
                  type = "perlin", seed = 1)
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "syn_noise")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "perlin")
})
