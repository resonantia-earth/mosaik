test_that("mdf_zonal summarises within zones", {

  # two patches: a 2x2 block of 1s and a single cell, on a 5x5 grid
  v <- c(1, 1, 0, 0, 0,
         1, 1, 0, 0, 0,
         0, 0, 0, 0, 0,
         0, 0, 0, 1, 0,
         0, 0, 0, 0, 0)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(fg = v))
  m <- mdf_componentise(m, connectivity = 8L, layer = "fg", add = "cc")
  # componentise labels the background as a component too, so it must be
  # masked out before it counts as a zone of its own
  m <- mdf_mask(m, by = "fg", background = 0, layer = "cc")

  # patch size onto every cell of its patch
  out <- mdf_zonal(m, by = "cc", fun = "n", layer = "cc", add = "size")
  sizes <- msk_pull(out, "size")

  expect_equal(sizes[1], 4)     # the 2x2 block
  expect_equal(sizes[19], 1)    # the lone cell
  expect_true(is.na(sizes[3]))  # background belongs to no zone
})


test_that("mdf_zonal 'any' detects a marker inside a zone", {

  # one 3-cell horizontal run; a marker layer flags only its middle cell
  v <- c(0, 0, 0, 0, 0,
         0, 1, 1, 1, 0,
         0, 0, 0, 0, 0,
         0, 1, 0, 0, 0,
         0, 0, 0, 0, 0)
  marker <- c(0, 0, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 0, 0,
              0, 0, 0, 0, 0,
              0, 0, 0, 0, 0)

  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(fg = v, marker = marker))
  m <- mdf_componentise(m, connectivity = 8L, layer = "fg", add = "cc")

  out <- mdf_zonal(m, by = "cc", fun = "any", layer = "marker", add = "hit")
  hits <- msk_pull(out, "hit")

  # every cell of the run carries the verdict, not just the marked one
  expect_equal(hits[7], 1)
  expect_equal(hits[8], 1)
  expect_equal(hits[9], 1)
  # the separate lone patch holds no marker
  expect_equal(hits[17], 0)
})


test_that("mdf_zonal counts distinct neighbouring labels", {

  # two blocks separated by a one-cell gap, bridged by the middle column
  #  A A . B B
  #  A A x B B      x = the connector cell, touching both blocks
  #  A A . B B
  blocks <- c(1, 1, 0, 1, 1,
              1, 1, 0, 1, 1,
              1, 1, 0, 1, 1)
  conn   <- c(0, 0, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 0, 0)

  m <- mosaik(extent = c(0, 5, 0, 3), res = 1,
              vals = list(blocks = blocks, conn = conn))
  m <- mdf_componentise(m, connectivity = 8L, layer = "blocks", add = "bid")
  m <- mdf_componentise(m, connectivity = 8L, layer = "conn", add = "cid")
  # background is a component as far as componentise is concerned; mask it off
  # both label layers so only real blocks and the real connector are zones
  m <- mdf_mask(m, by = "blocks", background = 0, layer = "bid")
  m <- mdf_mask(m, by = "conn", background = 0, layer = "cid")

  # the connector zone looks outward at the block IDs
  out <- mdf_zonal(m, by = "cid", fun = "n_distinct", neighbours = TRUE,
                   layer = "bid", add = "touch")
  touch <- msk_pull(out, "touch")

  # the connector touches two distinct blocks -> bridge
  expect_equal(touch[8], 2)
})


test_that("mdf_zonal neighbours excludes the zone's own cells", {

  # a single 3x3 block; its neighbourhood is all background (value 0)
  v <- c(0, 0, 0, 0, 0,
         0, 1, 1, 1, 0,
         0, 1, 1, 1, 0,
         0, 1, 1, 1, 0,
         0, 0, 0, 0, 0)
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(fg = v))
  m <- mdf_componentise(m, connectivity = 8L, layer = "fg", add = "cc")

  # summarising the fg layer over the neighbourhood must see only 0s,
  # never the block's own 1s
  out <- mdf_zonal(m, by = "cc", fun = "max", neighbours = TRUE,
                   layer = "fg", add = "nbr")
  expect_equal(msk_pull(out, "nbr")[7], 0)
})


test_that("mdf_zonal accepts a function and records a recipe step", {

  v <- c(1, 1, 0, 0,
         1, 1, 0, 0,
         0, 0, 2, 2,
         0, 0, 2, 2)
  m <- mosaik(extent = c(0, 4, 0, 4), res = 1, vals = list(z = v))

  # zones taken straight from the value layer (labels need not be components)
  out <- mdf_zonal(m, by = "z", fun = function(x) max(x) * 10,
                   layer = "z", add = "s")
  expect_equal(msk_pull(out, "s")[1], 10)
  expect_equal(msk_pull(out, "s")[11], 20)

  # recipe path: obj = NULL records instead of executing
  rec <- mdf_zonal(by = "z", fun = "n", layer = "z", add = "s")
  expect_s4_class(rec, "mosaik")
  expect_true(length(rec@provenance) == 1)
  expect_equal(names(rec@provenance)[1], "mdf_zonal")

  # and replaying it reproduces the direct call
  replayed <- mdf(m, rec)
  expect_equal(msk_pull(replayed, "s"), msk_pull(mdf_zonal(m, by = "z", fun = "n",
                                                           layer = "z", add = "s"), "s"))
})
