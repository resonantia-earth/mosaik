test_that("msk_select keeps specified layers", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(1, 25), b = rep(2, 25), c = rep(3, 25)))
  s <- msk_select(m, a, c)
  expect_equal(msk_names(s), c("a", "c"))
  expect_equal(msk_pull(s, "a"), rep(1, 25))
  expect_equal(msk_pull(s, "c"), rep(3, 25))
})

test_that("msk_select errors on missing layers", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(1, 25)))
  expect_error(msk_select(m, z), "none of the requested layers")
})

test_that("msk_add brings named layers in from a second mosaik", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(1, 25)))
  other <- mosaik(extent = c(0, 5, 0, 5), res = 1,
                  vals = list(b = rep(2, 25), c = rep(3, 25)))
  added <- msk_add(m, other, b)
  expect_equal(msk_names(added), c("a", "b"))
  expect_equal(msk_pull(added, "b"), rep(2, 25))
  # no layers named takes all of them
  expect_equal(msk_names(msk_add(m, other)), c("a", "b", "c"))
})

test_that("msk_add renames and refuses to overwrite", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(1, 25)))
  other <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(9, 25)))
  expect_error(msk_add(m, other, a), "already present")
  renamed <- msk_add(m, other, a, rename = "a2")
  expect_equal(msk_names(renamed), c("a", "a2"))
  expect_equal(msk_pull(renamed, "a"), rep(1, 25))
  expect_equal(msk_pull(renamed, "a2"), rep(9, 25))
})

test_that("msk_add rejects a mosaik on a different grid", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(1, 25)))
  wrongDims <- mosaik(extent = c(0, 5, 0, 5), res = 0.5,
                      vals = list(b = rep(2, 100)))
  wrongExtent <- mosaik(extent = c(0, 10, 0, 10), res = 2,
                        vals = list(b = rep(2, 25)))
  expect_error(msk_add(m, wrongDims, b), "same dimensions")
  expect_error(msk_add(m, wrongExtent, b), "same extent")
  expect_error(msk_add(m, m, z), "not found")
})

test_that("msk_add carries the added layer's categories", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(1, 25)))
  other <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(b = rep(2, 25)))
  other@categories$b <- list(role = "temperature")
  added <- msk_add(m, other, b)
  expect_equal(added@categories$b$role, "temperature")
})

test_that("mdf_filter masks cells", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = c(1:25)))
  f <- mdf_filter(m, cover > 20)
  vals <- msk_pull(f, "cover")
  expect_equal(sum(!is.na(vals)), 5)
  expect_true(all(vals[!is.na(vals)] > 20))
})

test_that("mdf_filter works with multi-layer predicates", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = 1:25, b = 25:1))
  f <- mdf_filter(m, a > b, layer = "a")
  vals <- msk_pull(f, "a")
  # cells where a <= b should be NA
  orig_a <- 1:25
  orig_b <- 25:1
  expect_equal(sum(!is.na(vals)), sum(orig_a > orig_b))
})

test_that("msk_pull finds layer values", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = 1:25))
  expect_equal(msk_pull(m, "cover"), 1:25)
})

test_that("msk_pull errors on missing layer", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1, 25)))
  expect_error(msk_pull(m, "nonexistent"), "not found")
})

test_that("derive computes class-level metric from mosaik", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(c(1L, 2L, 3L, 1L, 2L), 20)))
  m <- msr_area(m, scale = "class")
  m <- msr_perimeter(m, scale = "class")
  m <- msr(m, equation = "perimeter.class / area.class",
              label = "edge_density")
  expect_true("edge_density" %in% names(m@categories$cover))
  expect_equal(length(m@categories$cover$edge_density),
               length(m@categories$cover$gid))
})

test_that("derive computes mixed-scale metric", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(c(1L, 2L, 3L, 1L, 2L), 20)))
  m <- msr_area(m, scale = "class")
  m <- msr_area(m, scale = "landscape")
  m <- msr(m, equation = "area.class / area.landscape * 100",
              label = "prop_area")
  expect_true("prop_area" %in% names(m@categories$cover))
  expect_equal(sum(m@categories$cover$prop_area), 100)
})

test_that("derive errors on missing metric", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1L, 25)))
  expect_error(msr(m, equation = "area.class / area.landscape",
                      label = "x"), "not found")
})

test_that("derive errors on bad variable format", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1L, 25)))
  m <- msr_area(m, scale = "class")
  expect_error(msr(m, equation = "area", label = "x"),
               "metric.scale notation")
})

test_that("derive records provenance", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(c(1L, 2L), 50)))
  m <- msr_area(m, scale = "class")
  m <- msr_perimeter(m, scale = "class")
  m <- msr(m, equation = "perimeter.class / area.class",
              label = "edge_density")
  last_prov <- m@provenance[[length(m@provenance)]]
  expect_equal(names(last_prov)[1], "msr")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$equation, "perimeter.class / area.class")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$label, "edge_density")
})

test_that("derive with distance.cell computes GYRATE (mean centroid distance)", {
  # 10x10 grid, one 3x3 foreground patch
  vals <- rep(0L, 100)
  vals[c(23, 24, 25, 33, 34, 35, 43, 44, 45)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  m <- mdf_componentise(m)
  m <- mdf_centroid(m, add = "centroids")
  m <- mdf_distance(m, source = "centroids", add = "_distance")
  m <- msr(m, equation = "mean(distance.cell)", label = "gyrate")
  # result should be per-patch (2 patches: background + foreground)
  expect_true("gyrate" %in% names(m@patches))
  expect_equal(length(m@patches$gyrate), 2)
  # all values should be finite and non-negative
  expect_true(all(m@patches$gyrate >= 0))
})

test_that("derive with distance.cell computes max distance (CIRCLE)", {
  vals <- rep(0L, 100)
  vals[c(23, 24, 25, 33, 34, 35, 43, 44, 45)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  m <- mdf_componentise(m)
  m <- mdf_centroid(m, add = "centroids")
  m <- mdf_distance(m, source = "centroids", add = "_distance")
  m <- msr(m, equation = "max(distance.cell)", label = "circle")
  expect_true("circle" %in% names(m@patches))
  # max should be >= mean
  m <- msr(m, equation = "mean(distance.cell)", label = "gyrate")
  expect_true(all(m@patches$circle >= m@patches$gyrate))
})

test_that("derive .cell errors when internal layer missing", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(1L, 100)))
  expect_error(msr(m, equation = "mean(distance.cell)", label = "x"),
               "_distance.*not found")
})

test_that("msk_struct creates struct", {
  s <- msk_struct()
  expect_s4_class(s, "struct")
  expect_equal(dim(s@pattern), c(3, 3))
})

test_that("msk_struct custom pattern", {
  mat <- matrix(c(0,0,0,0,1,0,0,0,0), 3, 3)
  s <- msk_struct(custom = mat)
  expect_equal(s@pattern, mat)
})

test_that("msk_struct shapes differ and honour width and height", {
  expect_false(identical(msk_struct("disc", 5, 5)@pattern,
                         msk_struct("square", 5, 5)@pattern))
  expect_false(identical(msk_struct("disc", 5, 5)@pattern,
                         msk_struct("cross", 5, 5)@pattern))

  # width and height are respected, no recycling of a 3x3 pattern
  expect_equal(dim(msk_struct("disc", width = 7, height = 7)@pattern), c(7, 7))
  expect_equal(dim(msk_struct("rectangle", width = 7, height = 3)@pattern), c(3, 7))
  expect_silent(msk_struct("disc", width = 5, height = 5))
})

test_that("msk_struct disc follows the mmand definition", {
  # a cell counts when at least half of its area is inside the circle, so a
  # 3x3 disc is filled; the plus-shape is the 3x3 diamond
  expect_equal(msk_struct("disc", 3, 3)@pattern, matrix(1, nrow = 3, ncol = 3))

  d5 <- msk_struct("disc", 5, 5)@pattern
  expect_equal(d5, matrix(c(0, 1, 1, 1, 0,
                            1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1,
                            0, 1, 1, 1, 0), nrow = 5, ncol = 5, byrow = TRUE))
  expect_equal(d5[3, 3], 1)      # focal cell always set
  expect_equal(d5[1, 1], 0)      # corners cut

  expect_equal(msk_struct("disc", 7, 7)@pattern,
               matrix(c(0, 0, 1, 1, 1, 0, 0,
                        0, 1, 1, 1, 1, 1, 0,
                        1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1,
                        0, 1, 1, 1, 1, 1, 0,
                        0, 0, 1, 1, 1, 0, 0), nrow = 7, ncol = 7, byrow = TRUE))

  # discs grow monotonically with width, which granulometry depends on
  expect_true(sum(msk_struct("disc", 3, 3)@pattern) <
                sum(msk_struct("disc", 5, 5)@pattern))
  expect_true(sum(msk_struct("disc", 5, 5)@pattern) <
                sum(msk_struct("disc", 7, 7)@pattern))
})

test_that("msk_struct diamond is the manhattan ball", {
  # the historical default, relied on by mdf_dilate/mdf_erode/mdf_morph
  expect_equal(msk_struct("diamond", 3, 3)@pattern,
               matrix(c(0, 1, 0, 1, 1, 1, 0, 1, 0), nrow = 3, ncol = 3))

  expect_equal(msk_struct("diamond", 5, 5)@pattern,
               matrix(c(0, 0, 1, 0, 0,
                        0, 1, 1, 1, 0,
                        1, 1, 1, 1, 1,
                        0, 1, 1, 1, 0,
                        0, 0, 1, 0, 0), nrow = 5, ncol = 5, byrow = TRUE))
})

test_that("msk_struct stretches shapes to unequal width and height", {
  # the wider axis is spanned, rather than the largest ball fitting into both
  expect_equal(msk_struct("disc", width = 5, height = 3)@pattern,
               matrix(c(0, 1, 1, 1, 0,
                        1, 1, 1, 1, 1,
                        0, 1, 1, 1, 0), nrow = 3, ncol = 5, byrow = TRUE))
})

test_that("msk_struct square and cross are exact", {
  expect_equal(msk_struct("square", 5, 5)@pattern,
               matrix(1, nrow = 5, ncol = 5))
  expect_equal(msk_struct("rectangle", 5, 3)@pattern,
               matrix(1, nrow = 3, ncol = 5))
  expect_equal(msk_struct("cross", 5, 5)@pattern,
               matrix(c(0, 0, 1, 0, 0,
                        0, 0, 1, 0, 0,
                        1, 1, 1, 1, 1,
                        0, 0, 1, 0, 0,
                        0, 0, 1, 0, 0), nrow = 5, ncol = 5, byrow = TRUE))
})

test_that("msk_struct rejects even dimensions and unknown types", {
  expect_error(msk_struct(width = 4), "odd")
  expect_error(msk_struct(height = 4), "odd")
  expect_error(msk_struct(type = "blob"))
  # custom bypasses the shape checks entirely
  expect_silent(msk_struct(custom = matrix(1, nrow = 2, ncol = 4)))
})

test_that("msk_crop subsets extent and values", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = 1:100))
  cr <- mdf_crop(m, extent = c(2, 6, 3, 7))
  expect_equal(msk_dims(cr), c(4L, 4L))
  expect_equal(msk_extent(cr), c(2, 6, 3, 7))
  expect_equal(msk_ncells(cr), 16L)
})

test_that("msk_crop clamps to input extent", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(1L, 100)))
  cr <- mdf_crop(m, extent = c(-5, 15, -5, 15))
  expect_equal(msk_extent(cr), c(0, 10, 0, 10))
  expect_equal(msk_ncells(cr), 100L)
})

test_that("msk_crop errors on non-overlapping extent", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(1L, 100)))
  expect_error(mdf_crop(m, extent = c(20, 30, 20, 30)), "does not overlap")
})

test_that("msk_crop preserves multiple layers", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(a = 1:100, b = 101:200))
  cr <- mdf_crop(m, extent = c(0, 5, 0, 5))
  expect_equal(msk_names(cr), c("a", "b"))
  expect_equal(length(msk_pull(cr, "a")), msk_ncells(cr))
  expect_equal(length(msk_pull(cr, "b")), msk_ncells(cr))
})
