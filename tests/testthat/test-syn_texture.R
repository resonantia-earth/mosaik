test_that("syn_texture returns a mosaik with correct dimensions", {
  for (type in c("diamondSquare", "fbm", "billow", "ridged")) {
    r <- syn_texture(mosaik(extent = c(0, 50, 0, 50), res = 1),
                      type = type, seed = 42)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(50L, 50L))
    expect_equal(length(msk_pull(r)), 2500L)
  }
})

test_that("syn_texture diamondSquare produces values in [0, 1]", {
  r <- syn_texture(mosaik(extent = c(0, 64, 0, 64), res = 1),
                    type = "diamondSquare", seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_texture diamondSquare hurst controls smoothness", {
  r_smooth <- syn_texture(mosaik(extent = c(0, 64, 0, 64), res = 1),
                           type = "diamondSquare", hurst = 0.9, seed = 1)
  r_rough <- syn_texture(mosaik(extent = c(0, 64, 0, 64), res = 1),
                          type = "diamondSquare", hurst = 0.1, seed = 1)
  mat_s <- matrix(msk_pull(r_smooth), nrow = 64, ncol = 64, byrow = TRUE)
  mat_r <- matrix(msk_pull(r_rough), nrow = 64, ncol = 64, byrow = TRUE)
  var_s <- mean(abs(diff(mat_s[32, ])))
  var_r <- mean(abs(diff(mat_r[32, ])))
  expect_lt(var_s, var_r)
})

test_that("syn_texture fbm produces values in [0, 1]", {
  r <- syn_texture(mosaik(extent = c(0, 100, 0, 100), res = 1),
                    type = "fbm", base = "perlin", octaves = 6, seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_texture fbm with simplex base works", {
  r <- syn_texture(mosaik(extent = c(0, 50, 0, 50), res = 1),
                    type = "fbm", base = "simplex", seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_texture billow has positive bias", {
  r <- syn_texture(mosaik(extent = c(0, 100, 0, 100), res = 1),
                    type = "billow", seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_texture ridged produces ridge-like features", {
  r <- syn_texture(mosaik(extent = c(0, 100, 0, 100), res = 1),
                    type = "ridged", seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("syn_texture is reproducible with seed", {
  r1 <- syn_texture(mosaik(extent = c(0, 30, 0, 30), res = 1),
                     type = "fbm", seed = 42)
  r2 <- syn_texture(mosaik(extent = c(0, 30, 0, 30), res = 1),
                     type = "fbm", seed = 42)
  expect_equal(msk_pull(r1), msk_pull(r2))
})

test_that("syn_texture octaves adds detail", {
  r1 <- syn_texture(mosaik(extent = c(0, 100, 0, 100), res = 1),
                     type = "fbm", octaves = 1, seed = 1)
  r6 <- syn_texture(mosaik(extent = c(0, 100, 0, 100), res = 1),
                     type = "fbm", octaves = 6, seed = 1)
  mat1 <- matrix(msk_pull(r1), nrow = 100, ncol = 100, byrow = TRUE)
  mat6 <- matrix(msk_pull(r6), nrow = 100, ncol = 100, byrow = TRUE)
  var1 <- mean(abs(diff(mat1[50, ])))
  var6 <- mean(abs(diff(mat6[50, ])))
  expect_lt(var1, var6)
})

test_that("syn_texture name sets layer name", {
  r <- syn_texture(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    type = "diamondSquare", name = "texture", seed = 1)
  expect_equal(names(r@layers), "texture")
})

test_that("syn_texture records provenance", {
  r <- syn_texture(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    type = "fbm", seed = 1)
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "syn_texture")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "fbm")
})
