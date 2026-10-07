test_that("msr_area measures each class", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  a <- msr_area(obj = m)
  expect_true("area" %in% names(a@categories$cover))
  expect_equal(sum(a@categories$cover$area), 100)
})

test_that("on a layer of patch numbers, each class is a patch", {
  set.seed(42)
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  m <- mdf_componentise(m, add = "patch")
  a <- msr_area(m, layer = "patch")
  expect_equal(msk_table(a, "patch")$gid,
               sort(unique(msk_pull(m, "patch"))))
  expect_equal(sum(msk_table(a, "patch")$area), 100)
})

test_that("classes the map border cuts are measured and flagged by complete", {
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
    msr_area(layer = "patch") |>
    msr_perimeter(layer = "patch") |>
    msr_distance(routing = "straight", layer = "patch")
  pt <- msk_table(m, "patch")
  # no NA: the cut patch's values describe the part on the map
  expect_equal(pt$area, c(1, 2))
  expect_equal(pt$perimeter, c(2, 6))
  expect_equal(pt$distance["1", "2"], sqrt(8))
  expect_equal(pt$distance["2", "1"], sqrt(8))

  m <- msr(m, "complete.self", "whole", layer = "patch") |>
    msr("max(area.all[complete.all])", "largest", layer = "patch")
  expect_equal(msk_table(m, "patch")$whole, c(0, 1))
  expect_equal(msk_table(m, "patch")$largest, 2)
})

test_that("mdf_componentise numbers the patches of every class, 0 forms none", {
  v <- c(1, 1, 2,
         3, 1, 2,
         3, 3, 0)
  m <- mosaik(extent = c(0, 3, 0, 3), res = 1, vals = list(cover = v))
  m <- mdf_componentise(m, add = "patch")
  expect_true(is.na(msk_pull(m, "patch")[9]))
  expect_length(unique(na.omit(msk_pull(m, "patch"))), 3)
  # which class a patch came from is read through its cells
  m <- msr(m, "cover.self[1]", "source", layer = "patch")
  expect_equal(msk_table(m, "patch")$source, c(1, 2, 3))
})

test_that("msr_perimeter counts the edges to NA cells, not along the map border", {
  #   . . .
  #   . 1 .
  #   . . .
  m <- mosaik(extent = c(0, 3, 0, 3), res = 1,
              vals = list(cover = c(NA, NA, NA, NA, 1, NA, NA, NA, NA)))
  expect_equal(msk_table(msr_perimeter(m))$perimeter, 4)
  # the same cell at the border
  b <- mosaik(extent = c(0, 3, 0, 1), res = 1, vals = list(cover = c(1, 2, 2)))
  expect_equal(msk_table(msr_perimeter(b))$perimeter, c(1, 1))
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
  adj <- msk_table(msr_adjacency(obj = m))$adjacency

  # rook adjacencies, each pair counted from both sides:
  #   1-1: 4 undirected pairs -> 8;  2-2: 8 by symmetry
  #   1|2 seam: 2 undirected pairs -> 2 in each off-diagonal cell
  expect_equal(unname(adj["1", "1"]), 8)
  expect_equal(unname(adj["2", "2"]), 8)
  expect_equal(unname(adj["1", "2"]), 2)
  expect_equal(unname(adj["2", "1"]), 2)
})

test_that("the focus reads the adjacency matrix along the row of the class", {
  m <- landscape |> msr_adjacency(layer = "cover") |>
    msr("adjacency.self", "like", layer = "cover") |>
    msr("adjacency.self + sum(adjacency.others)", "paired", layer = "cover")
  cats <- msk_table(m, "cover")
  expect_equal(rownames(cats$adjacency), as.character(cats$gid))
  expect_equal(cats$like, unname(diag(cats$adjacency)))
  expect_equal(cats$paired, unname(rowSums(cats$adjacency)))
})

test_that("msr_adjacency counts the contacts between patches", {
  # forest(1) and grass(2) blocks side by side, water(3) below forest
  #   1 1 2 2
  #   1 1 2 2
  #   3 3 2 2
  v <- c(1, 1, 2, 2,
         1, 1, 2, 2,
         3, 3, 2, 2)
  m <- mosaik(extent = c(0, 4, 0, 3), res = 1, vals = list(cover = v))
  m <- mdf_componentise(m, add = "patch") |>
    msr_adjacency(connect = 4, layer = "patch")
  adj <- msk_table(m, "patch")$adjacency

  expect_equal(dim(adj), c(3, 3))
  expect_equal(adj, t(adj))
  expect_gt(adj["1", "2"], 0)
  expect_gt(adj["1", "3"], 0)
  expect_gt(adj["2", "3"], 0)
})

test_that("msr_distance measures between the classes of a layer", {
  # four corner patches of class 1
  vals <- c(1, 2, 2, 2, 1,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            1, 2, 2, 2, 1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_filter(m, cover == 1, add = "corner") |>
    mdf_replace(old = 0, new = NA, layer = "corner") |>
    mdf_componentise(layer = "corner", add = "patch") |>
    msr_distance(routing = "straight", layer = "patch")
  dmat <- msk_table(m, "patch")$distance
  expect_equal(dim(dmat), c(4, 4))
  expect_true(all(diag(dmat) == Inf))
  expect_equal(dmat, t(dmat))
  offdiag <- dmat[row(dmat) != col(dmat)]
  expect_true(all(offdiag > 0 & is.finite(offdiag)))

  # nearest neighbour of each corner: 4 cells away
  m <- msr(m, "min(distance.others)", "enn", layer = "patch")
  expect_equal(msk_table(m, "patch")$enn, rep(4, 4))

  # between land cover classes, in the order of the class table
  c2 <- msr_distance(m, routing = "straight", layer = "cover")
  expect_equal(rownames(msk_table(c2, "cover")$distance), c("1", "2"))
  expect_equal(msk_table(c2, "cover")$distance["1", "2"], 1)
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
  m <- mdf_filter(m, cover == 1, add = "one") |>
    mdf_replace(old = 0, new = NA, layer = "one") |>
    mdf_componentise(layer = "one", add = "patch")

  # the cheapest route crosses one cost-10 barrier cell plus cheap cells
  rs <- msr_distance(m, cost = "fric", routing = "cheapest", name = "effort",
                     layer = "patch")
  expect_equal(msk_table(rs, "patch")$effort[1, 2], 11)
  expect_null(msk_table(rs, "patch")$distance)

  # the name of a layer is not a free name
  expect_error(msr_distance(m, cost = "fric", name = "fric", layer = "patch"),
               "is a layer")
})

test_that("msr_distance returns NA across an impassable barrier", {
  # a full NA row seals the two class-1 blocks off from each other
  vals <- c(1,1,1,1,1,
            1,1,1,1,1,
            9,9,9,9,9,
            1,1,1,1,1,
            1,1,1,1,1)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = vals))
  m <- mdf_replace(m, old = c(1, 9), new = c(1, NA), add = "fric")
  m <- mdf_componentise(m, layer = "cover", add = "patch")

  r <- msr_distance(m, cost = "fric", routing = "cheapest", layer = "patch")
  mat <- msk_table(r, "patch")$distance
  expect_true(is.na(mat["1", "3"]))          # unreachable, distinct from a cost of 0
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

test_that("msr_dissimilarity weights each class's adjacencies", {
  vals <- c(1, 1, 2, 2,
            1, 1, 2, 2,
            3, 3, 3, 3,
            3, 3, 3, 3)
  m <- mosaik(extent = c(0, 4, 0, 4), res = 1, vals = list(cover = vals))

  cmat <- matrix(c(0, 0.5, 1.0,
                    0.5, 0, 0.8,
                    1.0, 0.8, 0), nrow = 3,
                 dimnames = list(c("1", "2", "3"), c("1", "2", "3")))

  m <- msr_dissimilarity(m, contrast = cmat)
  d <- m@categories$cover$dissimilarity
  expect_length(d, 3)
  expect_true(all(d >= 0))
  # the total weighted edge length of the layer
  m <- msr(m, "sum(dissimilarity.all) / 2", "total")
  expect_equal(msk_table(m)$total, sum(d) / 2)
})

test_that("msr builds a metric from an earlier metric", {
  m <- landscape |>
    msr_area(layer = "cover") |>
    msr(equation = "area.self / sum(area.all) * 100", label = "pland",
        layer = "cover") |>
    msr(equation = "-sum(pland.all / 100 * log(pland.all / 100))",
        label = "shannon", layer = "cover")
  a <- m@categories$cover$area
  p <- a / sum(a)
  expect_equal(msk_table(m, "cover")$shannon, -sum(p * log(p)))
  # a value stored for the layer is read with .all
  m <- msr(m, "shannon.all * 2", "double", layer = "cover")
  expect_equal(msk_table(m, "cover")$double, 2 * msk_table(m, "cover")$shannon)
  expect_error(msr(m, "shannon.self", "z", layer = "cover"), "shannon.all")
})

test_that("the focus words decide where the result is stored", {
  m <- msr_area(landscape, layer = "cover")
  m <- msr(m, "area.self * 2", "twice", layer = "cover") |>
    msr("max(area.all)", "most", layer = "cover")
  expect_equal(msk_table(m, "cover")$twice, 2 * msk_table(m, "cover")$area)
  expect_equal(msk_table(m, "cover")$most, max(msk_table(m, "cover")$area))
  # .others leaves the class in focus out
  m <- msr(m, "sum(area.others)", "rest", layer = "cover")
  cats <- msk_table(m, "cover")
  expect_equal(cats$rest, sum(cats$area) - cats$area)

  # a result that does not fit stops instead of being recycled
  expect_error(msr(m, "area.self + area.others", "z", layer = "cover"),
               "values for class .* Reduce")
  expect_error(msr(m, "area.all", "z", layer = "cover"),
               "values for the layer, where one is stored")
})

test_that("msr refuses a label that cannot be stored or read back", {
  m <- msr_area(landscape, layer = "cover")
  expect_error(msr(m, "area.self * 2", "area", layer = "cover"), "reserved name")
  expect_error(msr(m, "area.self * 2", "gid", layer = "cover"), "reserved name")
  expect_error(msr(m, "area.self * 2", "canopy", layer = "cover"), "is a layer")
  expect_error(msr(m, "area.self * 2", "a.b", layer = "cover"), "without '.'")
  m <- msr(m, "area.self * 2", "twice", layer = "cover")
  expect_error(msr(m, "area.self", "twice", layer = "cover"), "already stored")
  expect_error(msr(m, "area.class", "z", layer = "cover"), "<name>.<focus>")
  expect_error(msr(m, "perimeter.self", "z", layer = "cover"), "no value 'perimeter'")
  expect_error(msr(m, "area.all_nothere", "z", layer = "cover"), "not in 'obj'")
})

test_that("a layer name reads the cells of the classes in focus", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest") |>
    msr(equation = "mean(canopy.self)", label = "height", layer = "cover")
  cv <- msk_pull(m, "canopy"); lv <- msk_pull(m, "cover")
  cats <- msk_table(m, "cover")
  expect_equal(cats$height, vapply(cats$gid, function(k) mean(cv[lv == k]),
                                   numeric(1)))
  # a layer is read before a stored value of the same name
  m <- msr_area(m, layer = "forest")
  m@categories$forest$canopy <- c(-1, -1)
  m <- msr(m, "mean(canopy.self)", "h", layer = "forest")
  expect_false(any(msk_table(m, "forest")$h < 0))
})

test_that("x and y are the coordinates of the cell centres", {
  m <- mosaik(extent = c(10, 13, 0, 4), res = 1,
              vals = list(cover = c(1, 1, 2,
                                    1, 1, 2,
                                    2, 2, 2,
                                    2, 2, 2)))
  m <- msr(m, "mean(x.self)", "cx", layer = "cover") |>
    msr("mean(y.self)", "cy", layer = "cover")
  expect_equal(msk_table(m)$cx, c(11, mean(c(12.5, 12.5, 10.5, 11.5, 12.5,
                                                   10.5, 11.5, 12.5))))
  expect_equal(msk_table(m)$cy[1], 3)
})

test_that("another layer is read by its classes with .all, by cells otherwise", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest") |>
    mdf_replace(old = 0, new = NA, layer = "forest") |>
    mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
    msr_area(layer = "patch") |>
    msr_area(layer = "cover")
  m <- msr(m, "max(area.all) / sum(area.all_cover)", "lpi", layer = "patch")
  expect_equal(msk_table(m, "patch")$lpi,
               max(msk_table(m, "patch")$area) / msk_ncells(m))

  # a class value of cover, through the cells of each patch
  m@categories$cover$cost <- seq_along(m@categories$cover$gid)
  m <- msr(m, "mean(cost.self_cover)", "cost", layer = "patch")
  forest <- m@categories$cover$cost[m@categories$cover$gid == 47]
  expect_equal(unique(msk_table(m, "patch")$cost), forest)

  m <- msr_distance(m, routing = "straight", layer = "cover")
  expect_error(msr(m, "min(distance.self_cover)", "z", layer = "patch"),
               "only be read as a whole")
  expect_error(msr(m, "mean(canopy.self_cover)", "z", layer = "patch"),
               "no _layer suffix")
})

test_that("a layer name with underscores is read whole", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest_2020") |>
    msr_area(layer = "forest_2020") |>
    msr(equation = "sum(area.all_forest_2020) * 2", label = "double",
        layer = "cover")
  expect_equal(msk_table(m, "cover")$double, 6720)
})

test_that("msr reads pi and other base constants as constants", {
  m <- msr_area(landscape, layer = "cover") |>
    msr(equation = "2 * sqrt(area.self / pi)", label = "diameter",
        layer = "cover")
  p <- msk_table(m, "cover")
  expect_equal(p$diameter, 2 * sqrt(p$area / pi))
})

test_that("results are stored per layer and dropped when the layer is rewritten", {
  m <- landscape |>
    mdf_filter(cover == 47, add = "forest") |>
    msr_perimeter(layer = "cover") |>
    msr_perimeter(layer = "forest") |>
    msr("sum(perimeter.all)", "total", layer = "cover") |>
    msr("sum(perimeter.all)", "total", layer = "forest")
  expect_false(msk_table(m, "cover")$total == msk_table(m, "forest")$total)

  e <- mdf_erode(m, layer = "forest")
  expect_length(msk_table(e, "forest"), 0)
  expect_null(msk_table(e, "forest")$perimeter)
  expect_error(msr(e, "perimeter.self", "z", layer = "forest"), "no value")
})

test_that("msr_perimeter and msr_dissimilarity read a non-square grid the right way round", {
  # 2 rows x 4 cols; on a square grid a swap of rows and columns is invisible
  #   1 1 2 2
  #   1 1 2 2
  m <- mosaik(extent = c(0, 4, 0, 2), res = 1,
              vals = list(cover = c(1, 1, 2, 2, 1, 1, 2, 2)))
  # the two classes share a seam of two edges
  expect_equal(msk_table(msr_perimeter(m))$perimeter, c(2, 2))
  cmat <- matrix(c(0, 0.5, 0.5, 0), 2, dimnames = list(c("1", "2"), c("1", "2")))
  expect_equal(msk_table(msr_dissimilarity(m, contrast = cmat))$dissimilarity,
               c(1, 1))
  # an edge between left and right neighbours is as long as a cell is high
  r <- mosaik(extent = c(0, 4, 0, 6), res = c(1, 3),
              vals = list(cover = c(1, 1, 2, 2, 1, 1, 2, 2)))
  expect_equal(msk_table(msr_perimeter(r, unit = "map"))$perimeter, c(6, 6))
})
