test_that("msr_area at landscape scale", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  a <- msr_area(obj = m, scale = "landscape")
  expect_true("area" %in% names(a@global))
  expect_equal(a@global$area, 100)
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
  a <- msr_area(obj = m, scale = "patch")
  expect_true("area" %in% names(a@patches))
  expect_equal(sum(a@patches$area), 100)
})

test_that("msr_number at class scale counts classes", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:3, 100, replace = TRUE)))
  n <- msr_number(obj = m, scale = "class")
  expect_true("number" %in% names(n@global))
  expect_equal(n@global$number, 3L)
})

test_that("msr_number at patch scale counts patches per class", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:3, 100, replace = TRUE)))
  n <- msr_number(obj = m, scale = "patch")
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
  m <- msr_adjacency(m, scale = "patch", connect = 4)

  expect_true(all(c("adjacency", "regions") %in% names(m@patches)))
  adj <- m@patches$adjacency
  reg <- m@patches$regions

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
  m <- msr_adjacency(m, scale = "patch", connect = 8)

  reg <- m@patches$regions
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
  adj <- msr_adjacency(obj = m, type = "paired", count = "double")

  # double-counted rook adjacencies:
  #   1-1: 4 undirected pairs -> 8;  2-2: 8 by symmetry
  #   1|2 seam: 2 undirected pairs -> 2 in each off-diagonal cell
  expect_equal(unname(adj@global$adjacency["1", "1"]), 8)
  expect_equal(unname(adj@global$adjacency["2", "2"]), 8)
  expect_equal(unname(adj@global$adjacency["1", "2"]), 2)
  expect_equal(unname(adj@global$adjacency["2", "1"]), 2)
})

test_that("msr_cost computes pairwise patch distance matrices", {
  # 5x5 grid with two classes: class 1 in corners, class 2 elsewhere
  vals <- c(1, 2, 2, 2, 1,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            1, 2, 2, 2, 1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- msr_cost(m, routing = "straight")

  expect_true("distance" %in% names(m@patches))
  expect_true("1" %in% names(m@patches$distance))

  # class 1 has 4 patches (4 corners), so 4x4 matrix
  dmat <- m@patches$distance[["1"]]
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

test_that("msr_cost handles single-patch class", {
  # class 3 appears only once
  vals <- c(1, 1, 2, 2, 3,
            1, 1, 2, 2, 2,
            1, 2, 2, 2, 2,
            2, 2, 2, 2, 2)
  m <- mosaik(extent = c(0, 5, 0, 4), res = 1, vals = list(cover = vals))
  m <- msr_cost(m, routing = "straight")

  # class 3 has 1 patch: 1x1 matrix with Inf
  dmat <- m@patches$distance[["3"]]
  expect_equal(nrow(dmat), 1)
  expect_equal(dmat[1, 1], Inf)
})

test_that("msr_cost scale = 'cell' registers an edge-distance surface", {
  vals <- rep(0L, 100)
  vals[c(34, 35, 36, 44, 45, 46, 54, 55, 56)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  m <- msr_cost(m, scale = "cell")
  expect_true("_distance" %in% names(m@layers))
  dv <- msk_pull(m, "_distance")
  # centre of the foreground block should have largest edge distance
  expect_true(dv[45] > dv[34])
})

test_that("msr_cost with a cost surface accumulates sum and max", {
  # two class-1 blocks with an expensive band between them
  vals <- c(1,1,3,2,2,
            1,1,3,2,2,
            9,9,9,9,9,
            1,1,3,2,2,
            1,1,3,2,2)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_replace(m, old = c(1, 2, 3, 9), new = c(1, 1, 1, 10), add = "fric")

  # sum: the cheapest route crosses one cost-10 barrier cell plus cheap cells
  rs <- msr_cost(m, cost = "fric", routing = "cheapest", accumulate = "sum",
                 layer = "cover")
  expect_equal(rs@patches$fric[["1"]][1, 2], 11)

  # max: the worst cell on that same path is the barrier itself
  rm <- msr_cost(m, cost = "fric", routing = "cheapest", accumulate = "max",
                 layer = "cover")
  expect_equal(rm@patches$fric[["1"]][1, 2], 10)

  # the surface is named after what it measures, not the function
  expect_true("fric" %in% names(rs@patches))
  expect_false("distance" %in% names(rs@patches))

  # patch IDs are global: class 1 and the other classes share one ID space
  ids1 <- rownames(rs@patches$fric[["1"]])
  ids3 <- rownames(rs@patches$fric[["3"]])
  expect_length(intersect(ids1, ids3), 0)
})

test_that("msr_cost returns NA across an impassable barrier", {
  # a full NA row seals the two class-1 blocks off from each other
  vals <- c(1,1,1,1,1,
            1,1,1,1,1,
            9,9,9,9,9,
            1,1,1,1,1,
            1,1,1,1,1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_replace(m, old = c(1, 9), new = c(1, NA), add = "fric")

  r <- msr_cost(m, cost = "fric", routing = "cheapest", accumulate = "sum",
                layer = "cover")
  mat <- r@patches$fric[["1"]]
  expect_true(is.na(mat[1, 2]))          # unreachable, distinct from a cost of 0
  expect_true(all(is.infinite(diag(mat))))  # self stays Inf
})

test_that("msr_cost scale = 'cell' with a cost surface names the layer after it", {
  vals <- rep(0L, 100)
  vals[c(34, 35, 36, 44, 45, 46, 54, 55, 56)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  m <- mdf_replace(m, old = c(0, 1), new = c(2, 1), add = "friction")
  m <- msr_cost(m, scale = "cell", cost = "friction", layer = "cover")

  # the surface is named after what it measures, not after the function
  expect_true("_friction" %in% names(m@layers))
  expect_false("_distance" %in% names(m@layers))

  fv <- msk_pull(m, "_friction")
  # the source cells themselves cost nothing to reach
  expect_equal(fv[45], 0)
  # and cost accumulates away from them
  expect_true(fv[1] > 0)
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
  m <- msr_area(m, scale = "patch")
  m <- msr_cost(m, routing = "straight")
  m <- suppressWarnings(msr(m, equation = "min(distance.patch)", label = "enn"))

  # one ENN value per patch
  expect_equal(length(m@patches$enn), length(m@patches$patch))

  # class 1 has 4 corner patches, each nearest neighbour is 4 cells away
  cls1_idx <- which(m@patches$class == 1)
  expect_equal(m@patches$enn[cls1_idx], rep(4, 4))

  # class 2 has 1 patch, no neighbours → Inf
  cls2_idx <- which(m@patches$class == 2)
  expect_equal(m@patches$enn[cls2_idx], Inf)
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
  expect_true("dissimilarity" %in% names(m@global))
  expect_true(m@global$dissimilarity > 0)
})
