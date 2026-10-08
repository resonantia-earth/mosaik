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
  f <- mdf_filter(m, cover > 20, value = TRUE)
  vals <- msk_pull(f, "cover")
  expect_equal(sum(!is.na(vals)), 5)
  expect_true(all(vals[!is.na(vals)] > 20))
})

test_that("mdf_filter works with multi-layer predicates", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = 1:25, b = 25:1))
  f <- mdf_filter(m, a > b, value = TRUE, layer = "a")
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

test_that("msr computes a value per class from the class in focus", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(c(1L, 2L, 3L, 1L, 2L), 20)))
  m <- msr_area(m)
  m <- msr_perimeter(m)
  m <- msr(m, equation = "perimeter.self / area.self", label = "density")
  expect_true("density" %in% names(m@categories$cover))
  expect_equal(length(m@categories$cover$density),
               length(m@categories$cover$gid))
})

test_that("msr relates the class in focus to all classes", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(c(1L, 2L, 3L, 1L, 2L), 20)))
  m <- msr_area(m)
  m <- msr(m, equation = "area.self / sum(area.all) * 100", label = "share")
  expect_equal(sum(m@categories$cover$share), 100)
})

test_that("msr errors on a value that was not measured", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1L, 25)))
  expect_error(msr(m, equation = "area.self / sum(area.all)", label = "z"),
               "no value 'area'")
})

test_that("msr errors on bad variable format", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1L, 25)))
  m <- msr_area(m)
  expect_error(msr(m, equation = "area", label = "z"), "<name>.<focus>")
})

test_that("msr records provenance", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(c(1L, 2L), 50)))
  m <- msr_area(m)
  m <- msr_perimeter(m)
  m <- msr(m, equation = "perimeter.self / area.self", label = "density")
  last_prov <- m@provenance[[length(m@provenance)]]
  expect_equal(names(last_prov)[1], "msr")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$equation,
               "perimeter.self / area.self")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$label, "density")
})

test_that(".cell finds the cell a point lies in", {
  m <- mosaik(extent = c(0, 4, 0, 3), res = 1, vals = list(v = 1:12))
  # cells are numbered row by row from the top-left corner
  expect_equal(.cell(m, x = 0.5, y = 2.5), 1)
  expect_equal(.cell(m, x = 3.5, y = 0.5), 12)
  expect_equal(msk_pull(m, "v")[.cell(m, x = 1.5, y = 1.5)], 6)
  # points on the outer edge belong to the last column or row
  expect_equal(.cell(m, x = c(4, 0), y = c(0, 3)), c(12, 1))
  expect_true(is.na(.cell(m, x = 5, y = 1)))
})

test_that("a layer read through the cells of each patch gives GYRATE and CIRCLE", {
  # 10x10 grid, one 3x3 foreground patch
  vals <- rep(0L, 100)
  vals[c(23, 24, 25, 33, 34, 35, 43, 44, 45)] <- 1L
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1, vals = list(cover = vals))
  m <- mdf_componentise(m, add = "patch")
  m <- mdf_centroid(m, layer = "patch", add = "centroids")
  m <- mdf_distance(m, source = "centroids", layer = "patch", add = "dist")
  m <- msr(m, equation = "mean(dist.self)", label = "gyrate", layer = "patch")
  m <- msr(m, equation = "max(dist.self)", label = "circle", layer = "patch")
  p <- msk_table(m, "patch")
  # only the foreground patch is numbered
  expect_length(p$gyrate, 1)
  expect_true(p$gyrate >= 0)
  expect_true(p$circle >= p$gyrate)

  # the same from the cell coordinates, without a distance layer
  m <- msr(m, "mean(sqrt((x.self - mean(x.self))^2 + (y.self - mean(y.self))^2))",
           label = "gyr", layer = "patch")
  expect_gt(msk_table(m, "patch")$gyr, 0)
})

test_that("a name that is neither a layer nor a stored value stops", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = rep(1L, 100)))
  expect_error(msr(m, equation = "mean(dist.self)", label = "z"),
               "no value 'dist'")
})

test_that("the focus groups the cells", {
  #   1 1 2 2       canopy  10 20 30 30
  #   1 1 2 2               10 20 30 30
  v <- c(1, 1, 2, 2,
         1, 1, 2, 2)
  h <- c(10, 20, 30, 30,
         10, 20, 30, 30)
  m <- mosaik(extent = c(0, 4, 0, 2), res = 1,
              vals = list(cover = v, canopy = h))
  m <- msr(m, equation = "mean(canopy.self)", label = "height", layer = "cover")
  expect_equal(msk_table(m, "cover")$height, c(15, 30))
  m <- msr(m, equation = "mean(canopy.all)", label = "mean", layer = "cover")
  expect_equal(msk_table(m, "cover")$mean, 22.5)
  m <- msr(m, equation = "mean(canopy.others)", label = "around", layer = "cover")
  expect_equal(msk_table(m, "cover")$around, c(30, 15))

  # a cell value and a class value of the class in focus
  m <- msr_area(m, unit = "cells", layer = "cover")
  m <- msr(m, equation = "sum(canopy.self) / area.self", label = "avg",
           layer = "cover")
  expect_equal(msk_table(m, "cover")$avg, c(15, 30))
})

test_that("within a patch, a patch value is that patch's own value", {
  v <- c(0, 0, 0, 0, 0, 0,
         0, 1, 1, 0, 1, 0,
         0, 1, 1, 0, 1, 0,
         0, 0, 0, 0, 0, 0)
  h <- v * c(0, 0, 0, 0, 0, 0,
             0, 10, 20, 0, 40, 0,
             0, 10, 20, 0, 40, 0,
             0, 0, 0, 0, 0, 0)
  m <- mosaik(extent = c(0, 6, 0, 4), res = 1,
              vals = list(forest = v, canopy = h))
  m <- mdf_componentise(m, layer = "forest", add = "patch") |>
    msr_area(unit = "cells", layer = "patch") |>
    msr(equation = "sum(canopy.self) / area.self", label = "height",
        layer = "patch")
  expect_equal(msk_table(m, "patch")$height, c(15, 40))
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

test_that("msk_label writes labels and colours into the class table", {
  f <- mdf_filter(landscape, cover == 47, add = "forest")
  f <- msk_label(f, data.frame(id = 1, label = "forest"), layer = "forest")
  cats <- msk_table(f, "forest")
  # a class without a label so far is labelled with its code
  expect_equal(cats$val, c("0", "forest"))
  expect_null(cats$colour)
  f <- msk_label(f, data.frame(id = 0, label = "open", colour = "grey"),
                 layer = "forest")
  cats <- msk_table(f, "forest")
  expect_equal(cats$val, c("open", "forest"))
  expect_equal(cats$colour, c("grey", NA))
  # labels are not measured values
  expect_length(.measured_names(cats), 0)
  expect_error(msk_label(f, data.frame(id = 5, label = "x"), layer = "forest"),
               "do not occur")
  expect_error(msk_label(f, data.frame(id = 1), layer = "forest"))
  expect_error(msk_label(f, data.frame(id = 1, label = "a", size = 2),
                         layer = "forest"))
  expect_equal(names(utils::tail(f@provenance, 1)[[1]]), "msk_label")
})

test_that("changing the grid warns about values measured before", {
  m <- msr_area(landscape, layer = "cover") |>
    msr("max(area.all)", "most", layer = "cover")
  expect_warning(cr <- mdf_crop(m, c(0, 300, 0, 280)), "cover \\(area, most\\)")
  # the values are kept, as they were: the area of the whole map
  expect_equal(sum(msk_table(cr, "cover")$area), msk_ncells(m) * prod(msk_res(m)))
  expect_equal(msk_table(cr, "cover")$most, msk_table(m, "cover")$most)
  expect_warning(mdf_pad(m, width = 1L), "mdf_pad")
  expect_warning(mdf_resize(m, factor = 2), "mdf_resize")
  # labels alone are no reason to warn, and plotting a window does not warn
  expect_silent(mdf_crop(landscape, c(0, 300, 0, 280)))
  expect_silent(msk_vis(m, .layer("cover"), window = c(0, 30, 0, 28)))
})

test_that("msk_table is a list of the per-class and overall values by name", {
  m <- msr_area(landscape, layer = "cover") |>
    msr_distance(routing = "straight", layer = "cover") |>
    msr("max(area.all)", "most", layer = "cover")
  t <- msk_table(m, "cover")
  expect_s3_class(t, "msk_table")
  expect_true(is.list(t))
  expect_equal(t$area, m@categories$cover$area)
  expect_equal(t$most, max(t$area))
  out <- capture.output(print(t))
  expect_true(any(grepl("layer cover \\| 10 classes", out)))
  expect_true(any(grepl("matrices +distance \\(10 x 10\\)", out)))
  expect_true(any(grepl("overall +most = 83100", out)))
  expect_true(any(grepl("forest and hedgerows", out)))
  # the colours are not printed as text
  expect_false(any(grepl("#1f5f2e", out)))
  # beyond n classes, the rest are counted
  out <- capture.output(print(t, n = 3))
  expect_true(any(grepl("7 more classes", out)))
  expect_false(any(grepl("forest and hedgerows", out)))
  expect_output(print(msk_table(landscape, "canopy")), "0 classes")
  expect_error(msk_table(m, "nothere"), "not found")
})

test_that("a long argument in the history wraps, aligned with the arguments", {
  op <- options(width = 80)
  on.exit(options(op))
  m <- msr_area(landscape, layer = "cover") |>
    msr(equation = "-sum(area.all / sum(area.all) * log(area.all / sum(area.all)))",
        label = "shannon", layer = "cover") |>
    mdf_filter(cover == 1 | cover == 11 | cover == 21 | cover == 24 |
                 cover == 27 | cover == 31 | cover == 35 | cover == 41,
               add = "open")
  out <- capture.output(show(m))
  # msr shows its label, not its equation
  expect_true(any(grepl("msr .* label = \"shannon\"", out)))
  expect_false(any(grepl("equation", out)))

  out <- out[grep("mdf_filter", out):length(out)]
  expect_gt(length(out), 1)
  expect_true(all(nchar(out) <= 80))
  # continuation lines start in the argument column, never with an operator
  col <- regexpr("expr", out[1])
  expect_true(all(regexpr("\\S", out[-1]) >= col))
  expect_false(any(grepl("^\\s+(\\||[-+*/])", out[-1])))
})

test_that("legend labels can be switched off", {
  f <- msk_label(mdf_filter(landscape, cover == 47, add = "forest"),
                 data.frame(id = c(0, 1), label = c("open", "forest")),
                 layer = "forest")
  expect_true(.layer("forest")$labels)
  expect_false(.layer("forest", labels = FALSE)$labels)
  expect_silent(msk_vis(f, .layer("forest", labels = FALSE)))
  expect_error(.layer("forest", labels = NA))
})
