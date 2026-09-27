test_that("mdf_binarise works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = sample(1:5, 25, replace = TRUE)))
  b <- mdf_binarise(obj = m, thresh = 3)
  vals <- msk_pull(b, "cover")
  expect_true(all(vals %in% c(0L, 1L)))
})

test_that("mdf_offset works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1L, 25)))
  o <- mdf_offset(obj = m, value = 10)
  expect_equal(unique(msk_pull(o, "cover")), 11)
})

test_that("mdf_scale works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(2, 25)))
  s <- mdf_scale(obj = m, range = c(0, 1))
  vals <- msk_pull(s, "cover")
  # all same value scales to single point
  expect_length(unique(vals), 1)
})

test_that("mdf_replace works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = c(rep(1L, 13), rep(2L, 12))))
  r <- mdf_replace(obj = m, old = 1L, new = 99L)
  vals <- msk_pull(r, "cover")
  expect_false(1L %in% vals)
  expect_true(99L %in% vals)
})

test_that("mdf_permute works", {
  set.seed(1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = sample(1:3, 25, replace = TRUE)))
  # invert maps values 1,2,3 -> 2,1,0 (max - val)
  p <- mdf_permute(obj = m, type = "revert")
  # revert reverses the value mapping: 1->3, 2->2, 3->1
  orig <- msk_pull(m, "cover")
  perm <- msk_pull(p, "cover")
  expect_equal(length(perm), length(orig))
  expect_equal(sort(unique(perm)), sort(unique(orig)))
})

test_that("mdf_categorise works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = runif(25, 0, 10)))
  cc <- mdf_categorise(obj = m, n = 3)
  vals <- msk_pull(cc, "cover")
  expect_true(all(vals %in% 1:3))
})

test_that("mdf_componentise works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = c(1,0,1,0,1, 0,0,0,0,0, 1,0,1,0,1,
                                    0,0,0,0,0, 1,0,1,0,1)))
  b <- mdf_binarise(obj = m, match = 1L)
  co <- mdf_componentise(obj = b)
  vals <- msk_pull(co, "cover")
  # each isolated pixel should get its own component id
  expect_true(length(unique(vals[vals > 0])) > 1)
})

test_that("mdf_range works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = 1:25))
  r <- mdf_range(obj = m, lower = 10, upper = 20)
  vals <- msk_pull(r, "cover")
  expect_true(all(is.na(vals) | (vals >= 10 & vals <= 20)))
})

test_that("mdf_resize works", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(1, 100)))
  r <- mdf_resize(obj = m, factor = 2)
  expect_equal(msk_dims(r), c(20L, 20L))
  expect_equal(msk_ncells(r), 400L)
})

test_that("mdf_layerise works", {
  set.seed(1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = sample(1:3, 25, replace = TRUE)))
  result <- mdf_layerise(obj = m)
  expect_s4_class(result, "mosaik")
  expect_equal(length(result@layers), 3)
  for(nm in names(result@layers)){
    vals <- msk_pull(result, nm)
    # each layer has only one unique non-NA value
    expect_equal(length(unique(vals[!is.na(vals)])), 1)
  }
})

test_that("mdf_mask works", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = 1:25))
  m <- mdf_binarise(obj = m, thresh = 12, add = "mask")
  masked <- mdf_mask(obj = m, by = "mask", layer = "cover")
  vals <- msk_pull(masked, "cover")
  expect_true(any(is.na(vals)))
  expect_true(any(!is.na(vals)))
})

test_that("layers from a second mosaik blend after msk_add", {
  m1 <- mosaik(extent = c(0, 5, 0, 5), res = 1,
               vals = list(cover = rep(1, 25)))
  m2 <- mosaik(extent = c(0, 5, 0, 5), res = 1,
               vals = list(cover = rep(2, 25)))
  bl <- msk_add(m1, m2, cover, rename = "other") |>
    mdf_blend(layers = c("cover", "other"), fun = "+")
  expect_equal(unique(msk_pull(bl, "cover")), 3)
})

test_that("mdf_blend selects named layers and honours add", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(1, 25), b = rep(2, 25), c = rep(4, 25)))
  # only the named layers take part
  bl <- mdf_blend(obj = m, layers = c("a", "c"), fun = "+")
  expect_equal(unique(msk_pull(bl, "a")), 5)
  # defaults to every layer, written to a new one, leaving the inputs intact
  all <- mdf_blend(obj = m, fun = "+", add = "total")
  expect_equal(unique(msk_pull(all, "total")), 7)
  expect_equal(unique(msk_pull(all, "a")), 1)
})

test_that("mdf_blend weights each layer", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(2, 25), b = rep(4, 25)))
  bl <- mdf_blend(obj = m, fun = "+", weights = c(0.5, 0.25))
  expect_equal(unique(msk_pull(bl, "a")), 2)
})

test_that("mdf_blend rejects a bad layer request", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(1, 25), b = rep(2, 25)))
  expect_error(mdf_blend(obj = m, layers = c("a", "nope")), "not found")
  single <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(1, 25)))
  expect_error(mdf_blend(obj = single), "nothing to blend")
})

test_that("a second layer is named, never passed as a mosaik", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = 1:25))
  other <- mdf_binarise(obj = m, thresh = 12)

  # passing a mosaik is refused, and the error says what to do instead
  expect_error(mdf_mask(m, by = other), "msk_add")
  expect_error(mdf_zonal(m, by = other, fun = "n"), "msk_add")
  expect_error(mdf_layerise(m, by = other), "msk_add")

  # a name that is not a layer of obj is refused too
  expect_error(mdf_mask(m, by = "nope"), "not found")

  # the supported route: bring the layer in, then name it
  combined <- msk_add(m, other, cover, rename = "mask")
  expect_s4_class(mdf_mask(combined, by = "mask", layer = "cover"), "mosaik")
})

test_that("mdf_dilate and mdf_erode work", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  b <- mdf_binarise(obj = m, thresh = 3)
  d <- mdf_dilate(obj = b)
  e <- mdf_erode(obj = b)
  expect_true(all(msk_pull(d, "cover") %in% c(0, 1)))
  expect_true(all(msk_pull(e, "cover") %in% c(0, 1)))
  # dilate should have >= as many 1s as original
  expect_gte(sum(msk_pull(d, "cover") == 1), sum(msk_pull(b, "cover") == 1))
})

test_that("mdf_distance works", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  b <- mdf_binarise(obj = m, thresh = 3)
  d <- mdf_distance(obj = b)
  vals <- msk_pull(d, "cover")
  expect_true(min(vals) == 0)
  expect_true(max(vals) > 0)
})

test_that("mdf_blend with custom function", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(1, 25), b = rep(2, 25)))
  rd <- mdf_blend(obj = m, fun = sum)
  expect_equal(unique(msk_pull(rd, "a")), 3)
})

test_that("mdf_zonify burns a single polygon", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(0L, 100)))
  # square polygon covering x=2..8, y=2..8
  poly <- matrix(c(2,2, 8,2, 8,8, 2,8, 2,2), ncol = 2, byrow = TRUE)
  z <- mdf_zonify(m, geom = poly, value = 5L)
  vals <- msk_pull(z, "cover")
  # cells inside should be 5, outside should be NA (default background)
  expect_true(any(vals == 5L, na.rm = TRUE))
  expect_true(any(is.na(vals)))
})

test_that("mdf_zonify burns multiple polygons", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(0L, 100)))
  p1 <- matrix(c(0,0, 5,0, 5,5, 0,5, 0,0), ncol = 2, byrow = TRUE)
  p2 <- matrix(c(5,5, 10,5, 10,10, 5,10, 5,5), ncol = 2, byrow = TRUE)
  z <- mdf_zonify(m, geom = list(p1, p2), value = c(1L, 2L))
  vals <- msk_pull(z, "cover")
  expect_true(1L %in% vals)
  expect_true(2L %in% vals)
})

test_that("mdf_zonify with background = NULL preserves existing values", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(99L, 100)))
  poly <- matrix(c(3,3, 7,3, 7,7, 3,7, 3,3), ncol = 2, byrow = TRUE)
  z <- mdf_zonify(m, geom = poly, value = 1L, background = NULL)
  vals <- msk_pull(z, "cover")
  # cells outside polygon should still be 99
  expect_true(99L %in% vals)
  # cells inside should be 1
  expect_true(1L %in% vals)
})

test_that("mdf_zonify errors on unclosed polygon", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(0L, 100)))
  poly <- matrix(c(2,2, 8,2, 8,8, 2,8), ncol = 2, byrow = TRUE)
  expect_error(mdf_zonify(m, geom = poly, value = 1L), "not closed")
})

test_that("mdf_distance with coords snap = TRUE", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(0L, 100)))
  pts <- matrix(c(5, 5), ncol = 2)
  d <- mdf_distance(m, coords = pts, snap = TRUE)
  vals <- msk_pull(d)
  expect_equal(min(vals), 0)
  expect_true(max(vals) > 0)
})

test_that("mdf_distance with coords snap = FALSE", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(0L, 100)))
  pts <- matrix(c(5, 5), ncol = 2)
  d <- mdf_distance(m, coords = pts, snap = FALSE)
  vals <- msk_pull(d)
  # exact mode: min > 0 because point (5,5) doesn't sit on a centroid
  expect_true(min(vals) > 0)
})

test_that("mdf_distance source = 'background' measures edge distance", {
  # 10x10 grid, foreground block in centre
  vals <- rep(0L, 100)
  vals[c(34, 35, 36, 44, 45, 46, 54, 55, 56)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  d <- mdf_distance(m, source = "background")
  dv <- msk_pull(d)
  # background cells should still be 0 distance (they ARE background)
  expect_equal(dv[1], 0)
  # centre of the foreground block should have largest distance
  expect_true(dv[45] > dv[34])
  # foreground cells at edge should have distance ~1
  expect_true(dv[34] > 0)
})

test_that("mdf_distance source = layer name measures per-patch distance", {
  # create a binary layer with one patch, componentise, then centroid
  vals <- rep(0L, 100)
  vals[c(23, 24, 25, 33, 34, 35, 43, 44, 45)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  m <- mdf_componentise(m, add = "patches")
  m <- mdf_centroid(m, layer = "patches", add = "centroids")
  d <- mdf_distance(m, source = "centroids", layer = "patches", add = "dist")
  dv <- msk_pull(d, "dist")
  # componentise labels both background (patch 1) and foreground (patch 2),
  # so all cells get a distance value
  # centroid of foreground patch should have distance ~0
  expect_true(dv[34] < 0.1)
  # corner cells of the foreground patch should have larger distance
  expect_true(dv[23] > dv[34])
  # all cells should have a numeric distance (two patches cover everything)
  expect_false(any(is.na(dv)))
})

test_that("mdf_distance source = layer name errors when layer missing", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(1L, 100)))
  expect_error(mdf_distance(m, source = "nonexistent"), "not found")
})

test_that("mdf_pad adds cells", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = 1:25))
  p <- mdf_pad(m, width = 2, value = 0L)
  expect_equal(msk_dims(p), c(9L, 9L))
  expect_equal(msk_extent(p), c(-2, 7, -2, 7))
})

test_that("mdf_pad removes cells", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = 1:25))
  t <- mdf_pad(m, width = -1)
  expect_equal(msk_dims(t), c(3L, 3L))
})

test_that("mdf_pad single side", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1L, 25)))
  p <- mdf_pad(m, width = 3, sides = "right", value = 0L)
  expect_equal(msk_dims(p), c(8L, 5L))
  expect_equal(msk_extent(p)[2], 8)
})

test_that("mdf_rotate 180 reverses values", {
  m <- mosaik(extent = c(0, 4, 0, 3), res = 1,
              vals = list(cover = 1:12))
  r <- mdf_rotate(m, angle = 180)
  expect_equal(msk_pull(r), 12:1)
  expect_equal(msk_dims(r), msk_dims(m))
})

test_that("mdf_rotate 90 swaps dimensions", {
  m <- mosaik(extent = c(0, 4, 0, 3), res = 1,
              vals = list(cover = 1:12))
  r <- mdf_rotate(m, angle = 90)
  expect_equal(msk_dims(r), c(3L, 4L))
})

test_that("mdf_perturb adds noise", {
  set.seed(42)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(10, 25)))
  p <- mdf_perturb(m, sd = 1)
  vals <- msk_pull(p)
  expect_false(all(vals == 10))
  expect_true(abs(mean(vals) - 10) < 2)  # should be close to 10
})

test_that("mdf_perturb preserves NAs", {
  set.seed(1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = c(rep(10, 20), rep(NA, 5))))
  p <- mdf_perturb(m, sd = 1)
  expect_equal(sum(is.na(msk_pull(p))), 5)
})

test_that("mdf_interpolate smooths values", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = as.numeric(1:25)))
  sm <- mdf_interpolate(m)
  # smoothing should reduce the range
  expect_lte(diff(range(msk_pull(sm), na.rm = TRUE)),
             diff(range(msk_pull(m))))
})

test_that("mdf_tesselate partitions grid", {
  vals <- rep(NA, 100)
  vals[12] <- 1L; vals[55] <- 2L; vals[88] <- 3L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = vals))
  te <- mdf_tesselate(m)
  v <- msk_pull(te)
  expect_equal(sort(unique(v)), c(1L, 2L, 3L))
  expect_false(any(is.na(v)))
})

test_that("recipe records steps without an input mosaik", {
  recipe <- mdf_binarise(thresh = 3, layer = "cover", add = "fg") |>
    mdf_componentise(layer = "fg", add = "cc")
  steps <- Filter(function(e) isTRUE(e$step), recipe@provenance)
  expect_length(steps, 2)
  expect_equal(names(recipe@provenance)[1], "mdf_binarise")
  expect_equal(steps[[1]]$wasGeneratedBy$withArguments$add, "fg")
  expect_length(recipe@layers, 0)
})

test_that("mdf replays a recipe onto a real mosaik", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  recipe <- mdf_binarise(thresh = 3, layer = "cover", add = "fg") |>
    mdf_componentise(layer = "fg", add = "cc")
  out <- mdf(m, recipe)
  expect_true(all(c("cover", "fg", "cc") %in% msk_names(out)))
  # same result as running the steps directly
  direct <- mdf_componentise(mdf_binarise(m, thresh = 3, layer = "cover", add = "fg"),
                             layer = "fg", add = "cc")
  expect_equal(msk_pull(out, "cc"), msk_pull(direct, "cc"))
})

test_that("recipe step references an earlier layer by name via mdf_mask", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  recipe <- mdf_binarise(thresh = 3, layer = "cover", add = "fg") |>
    mdf_mask(layer = "cover", by = "fg", add = "masked")
  out <- mdf(m, recipe)
  expect_true("masked" %in% msk_names(out))
})

test_that("mdf errors on a recipe with no steps", {
  empty <- mosaik(extent = c(0, 1, 0, 1), res = 1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1L, 25)))
  expect_error(mdf(m, empty), "no recorded steps")
})

# --- fields another package attached, on in-place overwrite ----------------

# a layer carrying a field of the kind another package attaches to its
# category entry, which mosaik does not interpret
.tagged <- function() {
  g <- syn_gradient(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    name = "terrain")
  g@categories$terrain <- list(tag = "from another package")
  g
}
.tag_of <- function(obj, layer) obj@categories[[layer]]$tag

test_that("value-only mdf_* keep attached fields on overwrite", {
  g <- .tagged()
  for (op in list(function(x) mdf_scale(x, range = c(0, 100), layer = "terrain"),
                  function(x) mdf_perturb(x, layer = "terrain"),
                  function(x) mdf_offset(x, value = 5, layer = "terrain"))) {
    expect_equal(.tag_of(op(g), "terrain"), "from another package")
  }
})

test_that("kind-changing mdf_* drop attached fields on overwrite", {
  g <- .tagged()
  expect_null(.tag_of(mdf_binarise(g, thresh = 0.5, layer = "terrain"),
                      "terrain"))
  b <- mdf_binarise(g, thresh = 0.5, layer = "terrain")
  expect_null(.tag_of(mdf_componentise(b, layer = "terrain"), "terrain"))
})

test_that("writing to a new layer (add) leaves it without attached fields", {
  s <- mdf_scale(.tagged(), range = c(0, 1), layer = "terrain", add = "scaled")
  # the source keeps its field; the fresh layer has none to inherit
  expect_equal(.tag_of(s, "terrain"), "from another package")
  expect_null(.tag_of(s, "scaled"))
})

test_that("a categorical write starts a fresh category entry", {
  g <- .tagged()
  m <- msk_set(g, "terrain", rep(1:2, each = 50), gid = 1:2,
               val = c("low", "high"))
  expect_null(.tag_of(m, "terrain"))
  expect_equal(m@categories$terrain$val, c("low", "high"))
})
