test_that("mosaik(rast = ...) imports SpatRaster", {
  skip_if_not_installed("terra")

  r <- terra::rast(ncols = 10, nrows = 10,
                   xmin = 0, xmax = 10, ymin = 0, ymax = 10, nlyrs = 2)
  terra::values(r) <- cbind(sample(1:5, 100, replace = TRUE),
                            runif(100, 0, 1))
  names(r) <- c("cover", "ndvi")

  m <- mosaik(rast = r)
  expect_s4_class(m, "mosaik")
  expect_equal(msk_dims(m), c(10L, 10L))
  expect_equal(as.numeric(msk_extent(m)), c(0, 10, 0, 10))
  expect_equal(msk_pull(m, "cover"), terra::values(r[[1]])[, 1])
  expect_equal(msk_pull(m, "ndvi"), terra::values(r[[2]])[, 1], tolerance = 1e-10)
})

test_that("msk_terra roundtrip preserves data", {
  skip_if_not_installed("terra")

  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE),
                          ndvi = runif(100, 0, 1)))

  r <- suppressWarnings(msk_terra(m))
  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::nlyr(r), 2)
  expect_equal(names(r), c("cover", "ndvi"))

  m2 <- mosaik(rast = r)
  expect_equal(msk_pull(m2, "cover"), msk_pull(m, "cover"))
  expect_equal(msk_pull(m2, "ndvi"), msk_pull(m, "ndvi"), tolerance = 1e-10)
})

test_that("msk_terra warns about provenance loss", {
  skip_if_not_installed("terra")

  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(1, 25)))
  m2 <- mdf_offset(obj = m, value = 5)
  expect_warning(msk_terra(m2), "provenance")
})

test_that("mosaik(rast = ..., group = TRUE) creates categories", {
  skip_if_not_installed("terra")

  r <- terra::rast(ncols = 5, nrows = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5)
  terra::values(r) <- sample(1:3, 25, replace = TRUE)
  names(r) <- "cover"

  m <- mosaik(rast = r, group = TRUE)
  expect_true("cover" %in% names(m@categories))
  expect_true(!is.null(m@categories$cover$gid))
})
