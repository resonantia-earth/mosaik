test_that("msr_area at landscape scale", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  a <- msr_area(obj = m, scale = "landscape")
  expect_true("area" %in% names(msk_global(a)))
  expect_equal(msk_global(a)$area, 100)
})

test_that("msr_area at class scale", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  a <- msr_area(obj = m, scale = "class")
  expect_true("area" %in% names(a@categories$cover))
  expect_equal(sum(a@categories$cover$area), 100)
})

test_that("msr_area at patch scale", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  m <- mdf_pad(m, value = 0)
  m <- mdf_componentise(m, add = "patch")
  a <- msr_area(obj = m, scale = "patch")
  expect_true("area" %in% names(msk_patches(a)))
  expect_equal(sum(msk_patches(a)$area), 100)
})

test_that("patches the map border cuts have NA values", {
  #   1 0 0 0 0
  #   0 0 0 0 0
  #   0 0 1 1 0
  #   0 0 0 0 0
  v <- c(1, 0, 0, 0, 0,
         0, 0, 0, 0, 0,
         0, 0, 1, 1, 0,
         0, 0, 0, 0, 0)
  m <- mosaik(extent = c(0, 5, 0, 4), res = 1, vals = list(cover = v))
  m <- mdf_componentise(m, add = "patch") |>
    msr_area(scale = "patch") |>
    msr_perimeter(scale = "patch") |>
    msr_distance(routing = "straight") |>
    msr_adjacency(scale = "patch")
  pt <- msk_patches(m)
  expect_equal(pt$area, c(NA, 2))
  expect_equal(pt$perimeter, c(NA, 6))
  # the cut patch's own row is NA; the way to it from a whole patch is measured
  d <- pt$distance[["1"]]
  expect_true(all(is.na(d["1", ])))
  expect_equal(d["2", "1"], sqrt(8))
  expect_true(all(is.na(pt$adjacency["1", ])))
})

test_that("patch-level measures need patches numbered by mdf_componentise", {
  m <- mosaik(extent = c(0, 4, 0, 1), res = 1,
              vals = list(cover = c(1, 0, 1, 1)))
  expect_error(msr_area(m, scale = "patch"), "mdf_componentise")
  expect_error(msr_perimeter(m, scale = "patch"), "mdf_componentise")
  expect_error(msr_number(m, scale = "class"), "mdf_componentise")
  expect_error(msr_distance(m), "mdf_componentise")
  expect_error(msr_adjacency(m, scale = "patch"), "mdf_componentise")

  # the record says where the numbers are; 0 forms no patch
  m <- mdf_pad(m, value = 0)
  m <- mdf_componentise(m, connectivity = 8L, add = "patch")
  rec <- msk_patches(m)
  expect_equal(rec$ids, "patch")
  expect_equal(rec$class, c(1, 1))
  expect_equal(msr_area(m, scale = "patch")@patches$cover$area, c(1, 2))
})

test_that("mdf_componentise numbers the patches of every class", {
  v <- c(1, 1, 2,
         3, 1, 2,
         3, 3, 0)
  m <- mosaik(extent = c(0, 3, 0, 3), res = 1, vals = list(cover = v))
  m <- mdf_componentise(m, add = "patch")
  expect_equal(msk_patches(m)$class, c(1, 2, 3))
  expect_true(is.na(msk_pull(m, "patch")[9]))
  n <- msr_number(m, scale = "class")
  expect_equal(msk_categories(n)$number, c(0L, 1L, 1L, 1L))
})

test_that("the patch record goes when its layers are rewritten", {
  m <- mosaik(extent = c(0, 4, 0, 1), res = 1,
              vals = list(cover = c(1, 0, 1, 1)))
  m <- mdf_pad(m, value = 0)
  m <- mdf_componentise(m, add = "patch")
  # rewriting the source layer drops the record
  expect_null(msk_patches(mdf_replace(m, old = 0, new = 2), layer = "cover"))
  # rewriting the layer that holds the numbers drops it too
  expect_null(msk_patches(mdf_replace(m, old = 1, new = 5, layer = "patch"),
                          layer = "cover"))
  # overwriting the source layer with the numbers records them on it
  o <- mdf_componentise(m)
  expect_equal(msk_patches(o)$ids, "cover")
  # removing the number layer keeps the values but not the link
  r <- msk_remove(msr_area(m, scale = "patch"), patch)
  expect_equal(msk_patches(r)$area, c(1, 2))
  expect_error(msr_perimeter(r, scale = "patch"), "mdf_componentise")
})

test_that("msr_number at landscape scale counts classes", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:3, 100, replace = TRUE)))
  n <- msr_number(obj = m, scale = "landscape")
  expect_true("number" %in% names(msk_global(n)))
  expect_equal(msk_global(n)$number, 3L)
})

test_that("msr_number at class scale counts patches per class", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:3, 100, replace = TRUE)))
  m <- mdf_componentise(m, add = "patch")
  n <- msr_number(obj = m, scale = "class")
  expect_true("number" %in% names(n@categories$cover))
  expect_equal(length(n@categories$cover$number), 3L)
  expect_true(all(n@categories$cover$number >= 1L))
})

test_that("msr_perimeter works", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:3, 100, replace = TRUE)))
  p <- msr_perimeter(obj = m, scale = "class")
  expect_true("perimeter" %in% names(p@categories$cover))
})

test_that("msr_adjacency works", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:3, 100, replace = TRUE)))
  adj <- msr_adjacency(obj = m, type = "like")
  expect_true("likeAdj" %in% names(adj@categories$cover))
})

test_that("msr_adjacency at patch scale returns adjacency and region matrices", {
  # forest(1) and grass(2) blocks side by side, water(3) below forest
  #   1 1 2 2
  #   1 1 2 2
  #   3 3 2 2
  v <- c(1, 1, 2, 2,
         1, 1, 2, 2,
         3, 3, 2, 2)
  m <- mosaik(extent = c(0, 4, 0, 3), res = 1, vals = list(cover = v))
  m <- mdf_pad(m, value = 0)
  m <- mdf_componentise(m, add = "patch")
  m <- msr_adjacency(m, scale = "patch", connect = 4)

  expect_true(all(c("adjacency", "regions") %in% names(msk_patches(m))))
  adj <- msk_patches(m)$adjacency
  reg <- msk_patches(m)$regions

  # one patch per class -> 3 x 3
  expect_equal(dim(adj), c(3, 3))
  expect_equal(rownames(adj), c("1", "2", "3"))

  # adjacency is symmetric; regions need not be
  expect_equal(adj, t(adj))

  # forest(1) touches grass(2) and water(3); grass(2) touches water(3)
  expect_gt(adj["1", "2"], 0)
  expect_gt(adj["1", "3"], 0)
  expect_gt(adj["2", "3"], 0)

  # no loops in this layout -> every touching pair is a single region
  expect_equal(reg["1", "2"], 1)
  expect_equal(reg["1", "3"], 1)
})

test_that("msr_adjacency at patch scale detects a loop via the region count", {
  # core = patch 1 (2x2 top-left); fragment = patch 2 loops off it and touches
  # the core in TWO separate places (top-right and bottom), NA elsewhere.
  #   1 1 2 .
  #   1 1 . 2
  #   2 . . 2
  #   2 2 2 2
  v <- c(1, 1, 2, NA,
         1, 1, NA, 2,
         2, NA, NA, 2,
         2, 2, 2, 2)
  m <- mosaik(extent = c(0, 4, 0, 4), res = 1, vals = list(cover = v))
  m <- mdf_pad(m)
  m <- mdf_componentise(m, connectivity = 8L, add = "patch")
  m <- msr_adjacency(m, scale = "patch", connect = 8)

  reg <- msk_patches(m)$regions
  # read fragment -> core: two distinct contact places = loop signature
  expect_equal(reg["2", "1"], 2)
})

test_that("msr_adjacency is correct on a non-square grid", {
  # 2 rows x 4 cols; class 1 on the left half, class 2 on the right.
  # This distinguishes the (nrow, ncol) argument order: on a square grid a
  # swap is invisible, here it corrupts the counts.
  #   1 1 2 2
  #   1 1 2 2
  vals <- c(1, 1, 2, 2,
            1, 1, 2, 2)
  m <- mosaik(extent = c(0, 4, 0, 2), res = 1, vals = list(cover = vals))
  adj <- msk_categories(msr_adjacency(obj = m, type = "paired"))$adjacency

  # rook adjacencies, each pair counted from both sides:
  #   1-1: 4 undirected pairs -> 8;  2-2: 8 by symmetry
  #   1|2 seam: 2 undirected pairs -> 2 in each off-diagonal cell
  expect_equal(unname(adj["1", "1"]), 8)
  expect_equal(unname(adj["2", "2"]), 8)
  expect_equal(unname(adj["1", "2"]), 2)
  expect_equal(unname(adj["2", "1"]), 2)
})

test_that("paired adjacency is a class table: diagonal is like, row sums are pairedSum", {
  m <- landscape |>
    msr_adjacency(type = "paired", layer = "cover") |>
    msr_adjacency(type = "like", layer = "cover") |>
    msr_adjacency(type = "pairedSum", layer = "cover")
  cats <- msk_categories(m, "cover")
  expect_equal(rownames(cats$adjacency), as.character(cats$gid))
  expect_equal(unname(diag(cats$adjacency)), unname(cats$likeAdj))
  expect_equal(unname(rowSums(cats$adjacency)), unname(cats$pairedSum))
})

test_that("msr_distance computes pairwise patch distance matrices", {
  # 5x5 grid with two classes: class 1 in corners, class 2 elsewhere
  vals <- c(1, 2, 2, 2, 1,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            1, 2, 2, 2, 1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_pad(m, value = 0)
  m <- mdf_componentise(m, add = "patch")
  m <- msr_distance(m, routing = "straight")

  expect_true("distance" %in% names(msk_patches(m)))
  expect_true("1" %in% names(msk_patches(m)$distance))

  # class 1 has 4 patches (4 corners), so 4x4 matrix
  dmat <- msk_patches(m)$distance[["1"]]
  expect_equal(nrow(dmat), 4)
  expect_equal(ncol(dmat), 4)

  # diagonal is Inf
  expect_true(all(diag(dmat) == Inf))

  # matrix is symmetric
  expect_equal(dmat, t(dmat))

  # all off-diagonal distances are positive and finite
  offdiag <- dmat[row(dmat) != col(dmat)]
  expect_true(all(offdiag > 0))
  expect_true(all(is.finite(offdiag)))
})

test_that("msr_distance handles single-patch class", {
  # class 3 appears only once
  vals <- c(1, 1, 2, 2, 3,
            1, 1, 2, 2, 2,
            1, 2, 2, 2, 2,
            2, 2, 2, 2, 2)
  m <- mosaik(extent = c(0, 5, 0, 4), res = 1, vals = list(cover = vals))
  m <- mdf_pad(m, value = 0)
  m <- mdf_componentise(m, add = "patch")
  m <- msr_distance(m, routing = "straight")

  # class 3 has 1 patch: 1x1 matrix with Inf
  dmat <- msk_patches(m)$distance[["3"]]
  expect_equal(nrow(dmat), 1)
  expect_equal(dmat[1, 1], Inf)
})

test_that("msr_distance with a cost surface sums the costs along the path", {
  # two class-1 blocks with an expensive band between them
  vals <- c(1,1,3,2,2,
            1,1,3,2,2,
            9,9,9,9,9,
            1,1,3,2,2,
            1,1,3,2,2)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_pad(m, value = 9)
  m <- mdf_replace(m, old = c(1, 2, 3, 9), new = c(1, 1, 1, 10), add = "fric")
  m <- mdf_componentise(m, layer = "cover", add = "patch")

  # the cheapest route crosses one cost-10 barrier cell plus cheap cells
  rs <- msr_distance(m, cost = "fric", routing = "cheapest", layer = "cover")
  expect_equal(msk_patches(rs)$fric[["1"]][1, 2], 11)

  # the surface is named after what it measures, not the function
  expect_true("fric" %in% names(msk_patches(rs)))
  expect_false("distance" %in% names(msk_patches(rs)))

  # patch numbers are unique across classes: class 1 and the other classes
  # share one numbering
  ids1 <- rownames(msk_patches(rs)$fric[["1"]])
  ids3 <- rownames(msk_patches(rs)$fric[["3"]])
  expect_length(intersect(ids1, ids3), 0)
})

test_that("msr_distance returns NA across an impassable barrier", {
  # a full NA row seals the two class-1 blocks off from each other
  vals <- c(1,1,1,1,1,
            1,1,1,1,1,
            9,9,9,9,9,
            1,1,1,1,1,
            1,1,1,1,1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_pad(m, value = 9)
  m <- mdf_replace(m, old = c(1, 9), new = c(1, NA), add = "fric")
  m <- mdf_componentise(m, layer = "cover", add = "patch")

  r <- msr_distance(m, cost = "fric", routing = "cheapest", layer = "cover")
  mat <- msk_patches(r)$fric[["1"]]
  expect_true(is.na(mat[1, 2]))          # unreachable, distinct from a cost of 0
  expect_true(all(is.infinite(diag(mat))))  # self stays Inf
})

test_that("mdf_distance with a cost layer gives the cost of reaching each cell", {
  vals <- rep(0L, 100)
  vals[c(34, 35, 36, 44, 45, 46, 54, 55, 56)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  m <- mdf_replace(m, old = c(0, 1), new = c(2, 1), add = "friction")
  m <- mdf_distance(m, cost = "friction", layer = "cover", add = "reach")

  fv <- msk_pull(m, "reach")
  # the source cells themselves cost nothing to reach
  expect_equal(fv[45], 0)
  # and cost accumulates away from them
  expect_true(fv[1] > 0)
  expect_error(mdf_distance(m, source = "cover", cost = "friction"),
               "only with source")
})

test_that("msr_dissimilarity at class scale", {
  vals <- c(1, 1, 2, 2,
            1, 1, 2, 2,
            3, 3, 3, 3,
            3, 3, 3, 3)
  m <- mosaik(extent = c(0, 4, 0, 4), res = 1, vals = list(cover = vals))

  cmat <- matrix(c(0, 0.5, 1.0,
                    0.5, 0, 0.8,
                    1.0, 0.8, 0), nrow = 3,
                 dimnames = list(c("1", "2", "3"), c("1", "2", "3")))

  m <- msr_dissimilarity(m, contrast = cmat, scale = "class")
  expect_true("dissimilarity" %in% names(m@categories$cover))
  expect_equal(length(m@categories$cover$dissimilarity), 3)
  # all values should be non-negative
  expect_true(all(m@categories$cover$dissimilarity >= 0))
})

test_that("derive with distance.patch produces per-patch ENN", {
  vals <- c(1, 2, 2, 2, 1,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            1, 2, 2, 2, 1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_pad(m, value = 0)
  m <- mdf_componentise(m, add = "patch")
  m <- msr_distance(m, routing = "straight")
  m <- suppressWarnings(msr(m, equation = "min(distance.patch)", label = "enn"))

  # one ENN value per patch
  expect_equal(length(msk_patches(m)$enn), length(msk_patches(m)$patch))

  # class 1 has 4 corner patches, each nearest neighbour is 4 cells away
  cls1_idx <- which(msk_patches(m)$class == 1)
  expect_equal(msk_patches(m)$enn[cls1_idx], rep(4, 4))

  # class 2 has 1 patch, no neighbours → Inf
  cls2_idx <- which(msk_patches(m)$class == 2)
  expect_equal(msk_patches(m)$enn[cls2_idx], Inf)
})

test_that("msr_dissimilarity at landscape scale", {
  vals <- c(1, 1, 2, 2,
            1, 1, 2, 2,
            3, 3, 3, 3,
            3, 3, 3, 3)
  m <- mosaik(extent = c(0, 4, 0, 4), res = 1, vals = list(cover = vals))

  cmat <- matrix(c(0, 0.5, 1.0,
                    0.5, 0, 0.8,
                    1.0, 0.8, 0), nrow = 3,
                 dimnames = list(c("1", "2", "3"), c("1", "2", "3")))

  m <- msr_dissimilarity(m, contrast = cmat, scale = "landscape")
  expect_true("dissimilarity" %in% names(msk_global(m)))
  expect_true(msk_global(m)$dissimilarity > 0)
})

test_that("msr builds a metric from an earlier metric", {
  m <- landscape |>
    msr_area(scale = "class", layer = "cover") |>
    msr_area(scale = "landscape", layer = "cover") |>
    msr(equation = "area.class / area.landscape * 100", label = "pland",
        layer = "cover") |>
    msr(equation = "-sum(pland.class / 100 * log(pland.class / 100))",
        label = "shannon", layer = "cover")
  p <- m@categories$cover$area / msk_global(m)$area
  expect_equal(msk_global(m)$shannon, -sum(p * log(p)))
})

test_that("msr refuses a label that already exists at its level", {
  m <- landscape |>
    msr_area(scale = "class", layer = "cover") |>
    msr_area(scale = "landscape", layer = "cover")
  expect_error(msr(m, equation = "area.class * 2", label = "area",
                   layer = "cover"), "already exists at class level")
  expect_error(msr(m, equation = "area.class * 2", label = "gid",
                   layer = "cover"), "already exists at class level")
  expect_error(msr(m, equation = "area.landscape * 2", label = "area",
                   layer = "cover"), "already exists at landscape level")
})

test_that("results are stored per layer and do not overwrite each other", {
  m <- landscape |>
    mdf_pad() |>
    mdf_filter(cover == 47, add = "forest") |>
    msr_perimeter(scale = "landscape", layer = "cover") |>
    msr_perimeter(scale = "landscape", layer = "forest") |>
    mdf_componentise(layer = "cover", add = "cover_patch") |>
    mdf_componentise(layer = "forest", add = "forest_patch") |>
    msr_area(scale = "patch", layer = "cover") |>
    msr_area(scale = "patch", layer = "forest")
  expect_false(msk_global(m, "cover")$perimeter == msk_global(m, "forest")$perimeter)
  expect_equal(sum(msk_patches(m, "cover")$area), 3360)
  # only the forest (1) forms patches in a binary layer
  expect_equal(sum(msk_patches(m, "forest")$area),
               sum(msk_pull(m, "forest") == 1, na.rm = TRUE))
  expect_equal(unique(msk_patches(m, "forest")$class), 1)
})

test_that("rewriting a layer drops what was measured on it", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest") |>
    mdf_componentise(layer = "forest", add = "patch") |>
    msr_area(scale = "patch", layer = "forest") |>
    msr_area(scale = "landscape", layer = "forest") |>
    mdf_erode(layer = "forest")
  expect_null(msk_patches(m, "forest"))
  expect_null(msk_global(m, "forest"))
})

test_that("msr combines class values of two layers with the same classes", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest") |>
    mdf_erode(layer = "forest", add = "core") |>
    msr_area(scale = "class", layer = "forest") |>
    msr_area(scale = "class", layer = "core") |>
    msr(equation = "area.class_core / area.class_forest", label = "share",
        layer = "core")
  cats <- msk_categories(m, "core")
  expect_equal(cats$share,
               cats$area / msk_categories(m, "forest")$area)
})

test_that("msr refuses class or patch values of layers that do not correspond", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest") |>
    mdf_erode(layer = "forest", add = "core") |>
    msr_area(scale = "class", layer = "cover") |>
    msr_area(scale = "class", layer = "forest") |>
    mdf_componentise(layer = "forest", add = "forest_patch") |>
    mdf_componentise(layer = "core", add = "core_patch") |>
    msr_area(scale = "patch", layer = "forest") |>
    msr_area(scale = "patch", layer = "core")
  expect_error(msr(m, equation = "area.class_cover / area.class", label = "x",
                   layer = "forest"), "classes .* differ")
  expect_error(msr(m, equation = "area.patch_core / area.patch", label = "x",
                   layer = "forest"), "patches differ")
  expect_error(msr(m, equation = "area.class_nothere", label = "x",
                   layer = "forest"), "not in 'obj'")
})

test_that("a layer name with underscores is read whole", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest_2020") |>
    msr_area(scale = "landscape", layer = "forest_2020") |>
    msr(equation = "area.landscape_forest_2020 * 2", label = "double",
        layer = "cover")
  expect_equal(msk_global(m, "cover")$double, 6720)
})

test_that("the level in the label decides where the result is stored", {
  # a binary layer has one class, so one value fits the classes and the
  # landscape; the label says which
  f <- landscape |>
    mdf_filter(cover == 47, add = "forest") |>
    mdf_replace(old = 0, new = NA, layer = "forest") |>
    msr_area(scale = "class", layer = "forest") |>
    msr_area(scale = "landscape", layer = "cover")
  expect_length(msk_categories(f, "forest")$gid, 1)
  f <- f |>
    msr(equation = "area.class / area.landscape_cover", label = "pland.class",
        layer = "forest") |>
    msr(equation = "-sum(pland.class * log(pland.class))",
        label = "shdi.landscape", layer = "forest")
  expect_false(is.null(msk_categories(f, "forest")$pland))
  expect_false(is.null(msk_global(f, "forest")$shdi))
  expect_null(msk_categories(f, "forest")$shdi)

  # a result that does not fit the level stops instead of being recycled
  m <- landscape |>
    msr_area(scale = "class", layer = "cover") |>
    mdf_componentise(layer = "cover", add = "patch") |>
    msr_area(scale = "patch", layer = "cover")
  expect_error(msr(m, "area.patch / sum(area.class)", "x.class", layer = "cover"),
               "one value for each class")
  expect_error(msr(m, "area.class", "x.landscape", layer = "cover"),
               "takes one value")
  expect_error(msr(m, "area.class", "x.y.class", layer = "cover"), "without '.'")

  # without a level, the length decides, as before
  m <- msr(m, "area.class * 2", "twice", layer = "cover")
  expect_equal(msk_categories(m, "cover")$twice, 2 * msk_categories(m, "cover")$area)
})

test_that("msr reads pi and other base constants as constants", {
  m <- landscape |>
    mdf_componentise(layer = "cover", add = "patch") |>
    msr_area(scale = "patch", layer = "cover") |>
    msr(equation = "2 * sqrt(area.patch / pi)", label = "diameter.patch",
        layer = "cover")
  p <- msk_patches(m, "cover")
  expect_equal(p$diameter, 2 * sqrt(p$area / pi))
})

test_that("msr_perimeter and msr_dissimilarity read a non-square grid the right way round", {
  # 2 rows x 4 cols; on a square grid a swap of rows and columns is invisible
  #   1 1 2 2
  #   1 1 2 2
  m <- mosaik(extent = c(0, 4, 0, 2), res = 1,
              vals = list(cover = c(1, 1, 2, 2, 1, 1, 2, 2)))
  # the two classes share a seam of two edges
  expect_equal(msk_categories(msr_perimeter(m, scale = "class"))$perimeter, c(2, 2))
  cmat <- matrix(c(0, 0.5, 0.5, 0), 2, dimnames = list(c("1", "2"), c("1", "2")))
  expect_equal(msk_categories(msr_dissimilarity(m, contrast = cmat))$dissimilarity,
               c(1, 1))
  # an edge between left and right neighbours is as long as a cell is high
  r <- mosaik(extent = c(0, 4, 0, 6), res = c(1, 3),
              vals = list(cover = c(1, 1, 2, 2, 1, 1, 2, 2)))
  expect_equal(msk_categories(msr_perimeter(r, scale = "class", unit = "map"))$perimeter,
               c(6, 6))
})
