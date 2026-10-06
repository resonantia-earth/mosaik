grid10 <- function() mosaik(extent = c(0, 10, 0, 10), res = 1)

test_that("a polygon covers the cells whose centre lies inside it", {
  sq <- data.frame(x = c(2, 8, 8, 2), y = c(2, 2, 8, 8), id = 5)
  v <- msk_pull(msk_rasterise(grid10(), geom = sq))
  expect_equal(sum(v == 5, na.rm = TRUE), 36)
  expect_true(all(is.na(v[v != 5 | is.na(v)])))
  expect_equal(v[.cell(grid10(), x = 5.5, y = 5.5)], 5)
})

test_that("a closed ring and an open one give the same polygon", {
  open <- data.frame(x = c(2, 8, 8, 2), y = c(2, 2, 8, 8), id = 1)
  closed <- rbind(open, open[1, ])
  expect_equal(msk_pull(msk_rasterise(grid10(), geom = open)),
               msk_pull(msk_rasterise(grid10(), geom = closed)))
})

test_that("a ring inside another is a hole, and one inside that an island", {
  g <- data.frame(x = c(0, 10, 10, 0,  2, 8, 8, 2,  4, 6, 6, 4),
                  y = c(0, 0, 10, 10,  2, 2, 8, 8,  4, 4, 6, 6),
                  id = 1, part = rep(1:3, each = 4))
  m <- grid10()
  v <- msk_pull(msk_rasterise(m, geom = g))
  expect_equal(v[.cell(m, x = 0.5, y = 0.5)], 1)   # outer ring
  expect_true(is.na(v[.cell(m, x = 3, y = 3)]))    # hole
  expect_equal(v[.cell(m, x = 5, y = 5)], 1)       # island
  expect_equal(sum(!is.na(v)), 100 - 36 + 4)
})

test_that("later geometries win where they overlap", {
  g <- data.frame(x = c(0, 6, 6, 0,  4, 10, 10, 4),
                  y = c(0, 0, 6, 6,  4, 4, 10, 10),
                  id = rep(c(1, 2), each = 4))
  m <- grid10()
  v <- msk_pull(msk_rasterise(m, geom = g))
  expect_equal(v[.cell(m, x = 5, y = 5)], 2)
})

test_that("points cover their cells and points outside are dropped", {
  p <- data.frame(x = c(0.5, 9.5, 20), y = c(0.5, 9.5, 20), id = c(1, 2, 3))
  m <- grid10()
  v <- msk_pull(msk_rasterise(m, geom = p, type = "point"))
  expect_equal(sum(!is.na(v)), 2)
  expect_equal(v[.cell(m, x = 9.5, y = 9.5)], 2)
})

test_that("a line covers the cells it passes through, part by part", {
  l <- data.frame(x = c(0.5, 9.5, 0.5, 9.5), y = c(0.5, 0.5, 5.5, 5.5),
                  id = 1, part = c(1, 1, 2, 2))
  m <- grid10()
  v <- msk_pull(msk_rasterise(m, geom = l, type = "line"))
  # two horizontal rows, and nothing in between: the parts are not joined
  expect_equal(sum(!is.na(v)), 20)
  expect_true(is.na(v[.cell(m, x = 5, y = 3)]))
})

test_that("rasterised points feed mdf_distance", {
  p <- data.frame(x = 5.5, y = 5.5, id = 1)
  m <- grid10() |>
    msk_rasterise(geom = p, type = "point", name = "pt") |>
    mdf_filter(pt == 1, add = "src") |>
    mdf_distance(layer = "src", add = "d")
  v <- msk_pull(m, "d")
  expect_equal(v[.cell(m, x = 5.5, y = 5.5)], 0)
  expect_equal(v[.cell(m, x = 8.5, y = 5.5)], 3)
})

test_that("malformed geometries are refused", {
  expect_error(msk_rasterise(grid10(), geom = data.frame(x = 1, y = 1)))
  expect_error(msk_rasterise(grid10(),
                             geom = data.frame(x = c(1, 2), y = c(1, 2), id = 1)),
               "three vertices")
  expect_error(msk_rasterise(grid10(), geom = data.frame(x = 1, y = 1, id = 1),
                             type = "line"), "two vertices")
})

test_that("msk_spaghettify numbers features and their parts", {
  skip_if_not_installed("sf")
  outer <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0), ncol = 2, byrow = TRUE)
  hole <- matrix(c(2, 2, 8, 2, 8, 8, 2, 8, 2, 2), ncol = 2, byrow = TRUE)
  small <- matrix(c(20, 20, 21, 20, 21, 21, 20, 20), ncol = 2, byrow = TRUE)
  x <- sf::st_sfc(sf::st_multipolygon(list(list(outer, hole), list(small))),
                  sf::st_polygon(list(outer)))
  s <- msk_spaghettify(x)
  expect_named(s, c("x", "y", "id", "part"))
  expect_equal(unique(s$id), c(1, 2))
  expect_equal(unique(s$part[s$id == 1]), c(1, 2, 3))
  expect_equal(unique(s$part[s$id == 2]), 1)

  # the hole survives the round trip into the grid
  m <- msk_rasterise(grid10(),
                     geom = s[s$id == 1, ])
  v <- msk_pull(m)
  expect_true(is.na(v[.cell(m, x = 5, y = 5)]))

  l <- sf::st_sfc(sf::st_linestring(matrix(c(0, 0, 5, 5, 9, 2), ncol = 2,
                                           byrow = TRUE)))
  expect_equal(unique(msk_spaghettify(l)$part), 1)
  expect_error(msk_spaghettify(c(x, l)), "one kind")
})
