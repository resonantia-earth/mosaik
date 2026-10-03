test_that("mosaik() creates valid object", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(1L, 100)))
  expect_s4_class(m, "mosaik")
  expect_equal(msk_extent(m), c(0, 10, 0, 10))
  expect_equal(msk_dims(m), c(10L, 10L))
  expect_equal(msk_res(m), c(1, 1))
  expect_equal(msk_ncells(m), 100L)
  expect_equal(msk_crs(m), NA_character_)
  expect_equal(msk_names(m), "cover")
  expect_equal(msk_pull(m, "cover"), rep(1L, 100))
})

test_that("mosaik() with CRS", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, crs = "+proj=longlat",
              vals = list(a = rep(0, 100)))
  expect_equal(msk_crs(m), "+proj=longlat")
})

test_that("mosaik() validates bad inputs", {
  expect_error(mosaik(extent = c(10, 0, 0, 10), res = 1,
                      vals = list(a = rep(1, 100))))
  expect_error(mosaik(extent = c(0, 10, 0, 10), res = 1,
                      vals = list(a = rep(1, 50))))
})

test_that("msk_pull default is first layer", {
  m <- mosaik(extent = c(0, 4, 0, 4), res = 1,
              vals = list(a = 1:16, b = 16:1))
  expect_equal(msk_pull(m), 1:16)
})

test_that("mosaik() constructor adds provenance", {
  m <- mosaik(extent = c(0, 2, 0, 2), res = 1,
              vals = list(x = rep(1, 4)))
  expect_length(msk_provenance(m), 1)
})

test_that("the class validity rejects inconsistent objects", {
  expect_s4_class(methods::new("mosaik", extent = c(0, 5, 0, 5), dims = c(5L, 5L),
                               layers = list(a = rep(1, 25)),
                               crs = NA_character_), "mosaik")
  expect_error(methods::new("mosaik", extent = c(0, 5, 0, 5), dims = c(5L, 5L),
                            layers = list(a = rep(1, 24)),
                            crs = NA_character_), "24 values")
  expect_error(methods::new("mosaik", extent = c(0, 5, 0, 5), dims = c(5L, 5L),
                            layers = list(a = rep(1, 25)),
                            categories = list(b = list(gid = 1)),
                            crs = NA_character_), "does not match any layer")
})

test_that("show method runs without error", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = sample(1:3, 25, replace = TRUE)))
  expect_output(show(m))
})
