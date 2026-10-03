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

  # all values are replaced at once: swapping 1 and 2 swaps them
  s <- msk_pull(mdf_replace(obj = m, old = c(1, 2), new = c(2, 1)), "cover")
  expect_equal(s, c(rep(2, 13), rep(1, 12)))
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

test_that("mdf_categorise by shares hits the shares, lowest values first", {
  set.seed(1)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(v = c(NA, runif(99))))
  cc <- mdf_categorise(m, shares = c(0.5, 0.3, 0.2), add = "cl")
  cl <- msk_pull(cc, "cl")
  v <- msk_pull(cc, "v")
  # 99 valued cells: round(49.5) = 50, round(79.2) = 79, 99
  expect_equal(tabulate(cl), c(50, 29, 20))
  expect_true(is.na(cl[1]))
  expect_true(max(v[cl %in% 1]) < min(v[cl %in% 2]))

  # shares that do not sum to 1 are refused, and so is more than one mode
  expect_error(mdf_categorise(m, shares = c(0.5, 0.3)), "sum to 1")
  expect_error(mdf_categorise(m, n = 2, shares = c(0.5, 0.5)), "exactly one")

  # equal values share a category, so a binary layer cannot be split 30/70
  b <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(v = rep(c(0, 1), each = 50)))
  expect_error(mdf_categorise(b, shares = c(0.3, 0.7)), "cannot be met")
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

test_that("a second layer is named, never passed as a mosaik", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = 1:25))
  other <- mdf_binarise(obj = m, thresh = 12)

  # passing a mosaik is refused, and the error says what to do instead
  expect_error(mdf_summarise(m, by = other, fun = "n"), "msk_add")

  # a name that is not a layer of obj is refused too
  expect_error(mdf_summarise(m, by = "nope", fun = "n"), "not found")

  # the supported route: bring the layer in, then name it
  combined <- msk_add(m, other, cover, rename = "zones")
  expect_s4_class(mdf_summarise(combined, by = "zones", fun = "n",
                                layer = "cover"), "mosaik")
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
  # centroid of the patch should have distance ~0
  expect_true(dv[34] < 0.1)
  # corner cells of the patch should have larger distance
  expect_true(dv[23] > dv[34])
  # only the patch is componentised, so only its cells get a distance
  expect_false(any(is.na(dv[vals == 1])))
  expect_true(all(is.na(dv[vals == 0])))
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
  steps <- recipe@provenance
  expect_length(steps, 2)
  expect_equal(names(steps[[1]]), "mdf_binarise")
  expect_equal(steps[[1]][[1]]$wasGeneratedBy$withArguments$add, "fg")
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

test_that("recipe step references an earlier layer by name", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  recipe <- mdf_binarise(thresh = 3, layer = "cover", add = "fg") |>
    mdf_summarise(by = "fg", fun = "n", layer = "cover", add = "masked")
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
  g <- drw_gradient(mosaik(extent = c(0, 10, 0, 10), res = 1),
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
  m <- mdf_categorise(g, n = 2, layer = "terrain")
  expect_null(.tag_of(m, "terrain"))
  expect_length(m@categories$terrain$val, 2)
})

# --- measurements inside a recipe ------------------------------------------

test_that("msr_* and msr() record into a recipe and replay with mdf()", {
  rec <- mdf_binarise(match = 47, layer = "cover", add = "forest") |>
    msr_area(scale = "class", layer = "forest") |>
    mdf_erode(layer = "forest", add = "core") |>
    msr_area(scale = "class", layer = "core") |>
    msr_area(scale = "landscape", layer = "core") |>
    msr(equation = "area.class / area.landscape * 100", label = "pcore",
        layer = "core")

  # nothing was computed while recording
  expect_length(rec@layers, 0)

  direct <- landscape |>
    mdf_binarise(match = 47, layer = "cover", add = "forest") |>
    msr_area(scale = "class", layer = "forest") |>
    mdf_erode(layer = "forest", add = "core") |>
    msr_area(scale = "class", layer = "core") |>
    msr_area(scale = "landscape", layer = "core") |>
    msr(equation = "area.class / area.landscape * 100", label = "pcore",
        layer = "core")

  replayed <- mdf(landscape, rec)
  expect_equal(replayed@categories, direct@categories)
  expect_equal(replayed@global, direct@global)
})

test_that("an argument recorded into a recipe is stored by value", {
  cls <- as.character(sort(unique(msk_pull(landscape, "cover"))))
  ct <- matrix(1, length(cls), length(cls), dimnames = list(cls, cls))
  diag(ct) <- 0
  rec <- msr_dissimilarity(contrast = ct, layer = "cover")
  expected <- msr_dissimilarity(landscape, contrast = ct, layer = "cover")
  rm(ct)
  expect_equal(mdf(landscape, rec)@categories$cover$dissimilarity,
               expected@categories$cover$dissimilarity)
})

test_that("mdf_morph with identity and max/min is dilation/erosion", {
  f <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest")
  s <- msk_struct("square", width = 3, height = 3)
  morph <- function(merge) {
    msk_pull(mdf_morph(f, struct = s, blend = "identity", merge = merge,
                       layer = "forest", add = "out"), "out")
  }
  expect_equal(morph("max"),
               msk_pull(mdf_dilate(f, struct = s, layer = "forest", add = "out"), "out"))
  expect_equal(morph("min"),
               msk_pull(mdf_erode(f, struct = s, layer = "forest", add = "out"), "out"))
})

test_that("mdf_morph accepts a kernel wider than it is tall", {
  s <- msk_struct("square", width = 5, height = 3)
  m <- mdf_morph(landscape, struct = s, blend = "identity", merge = "mean",
                 layer = "canopy", add = "smooth")
  expect_true("smooth" %in% msk_names(m))
})

test_that("mdf_morph refuses to rotate a kernel that is not square", {
  s <- msk_struct("square", width = 5, height = 3)
  expect_error(mdf_morph(landscape, struct = s, blend = "equal", merge = "all",
                         rotate = TRUE, layer = "canopy"),
               "needs a square")
})

test_that("mdf_morph respects the shape of the structuring element", {
  f <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest")
  d <- msk_struct("disc", width = 5, height = 5)
  expect_equal(
    msk_pull(mdf_morph(f, struct = d, blend = "identity", merge = "max",
                       layer = "forest", add = "out"), "out"),
    msk_pull(mdf_dilate(f, struct = d, layer = "forest", add = "out"), "out"))
})

test_that("mdf_match gives no result where the pattern does not fit the grid", {
  # a single forest cell in the corner and one in the middle; the corner one
  # cannot be tested with a 3 x 3 pattern and must not keep its input value
  v <- rep(0, 25); v[c(1, 13)] <- 1
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(f = v))
  iso <- msk_struct(custom = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), nrow = 3))
  r <- msk_pull(mdf_match(m, struct = iso, layer = "f", add = "r"), "r")
  expect_equal(which(!is.na(r)), 13)
  expect_true(all(is.na(r[c(1:5, 21:25)])))
})

test_that("an NA neighbour does not satisfy a pattern cell in mdf_match", {
  # an isolated 1 whose neighbours are 0, except one that is NA
  v <- rep(0, 25); v[13] <- 1; v[12] <- NA
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(f = v))
  iso <- msk_struct(custom = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), nrow = 3))
  r <- msk_pull(mdf_match(m, struct = iso, layer = "f", add = "r"), "r")
  expect_true(all(is.na(r)))
})

test_that("mdf_categorise takes breaks beyond the values", {
  m <- mosaik(extent = c(0, 4, 0, 1), res = 1, vals = list(v = c(0, 1.5, 2.5, 3)))
  cc <- msk_pull(mdf_categorise(m, breaks = 1:8, layer = "v"), "v")
  # [0, 1), [1, 2), [2, 3), [3, 4): one category each
  expect_equal(cc, c(1, 2, 3, 4))
})
