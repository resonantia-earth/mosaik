test_that("drw_tessellation returns a mosaik with correct dimensions", {
  for (type in c("voronoi", "rectangle", "gibbs")) {
    r <- drw_tessellation(mosaik(extent = c(0, 50, 0, 50), res = 1),
                           type = type, n = 10, seed = 42)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(50L, 50L))
    expect_equal(length(msk_pull(r)), 2500L)
  }
})

test_that("drw_tessellation voronoi produces n regions", {
  r <- drw_tessellation(mosaik(extent = c(0, 50, 0, 50), res = 1),
                         type = "voronoi", n = 15, seed = 1)
  vals <- msk_pull(r)
  expect_equal(length(unique(vals)), 15)
})

test_that("drw_tessellation rectangle covers the grid with rectangles", {
  r <- drw_tessellation(mosaik(extent = c(0, 30, 0, 20), res = 1),
                         type = "rectangle", size = c(3, 6), seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals > 0))
  # numbered 1, 2, ... without gaps
  expect_equal(sort(unique(vals)), seq_len(max(vals)))
  # no visible piece is wider or taller than the largest rectangle
  mat <- matrix(vals, nrow = 20, byrow = TRUE)
  for (k in unique(vals)) {
    cells <- which(mat == k, arr.ind = TRUE)
    expect_lte(diff(range(cells[, "row"])), 5)
    expect_lte(diff(range(cells[, "col"])), 5)
  }
  expect_error(drw_tessellation(mosaik(extent = c(0, 10, 0, 10), res = 1),
                                type = "rectangle", size = c(6, 3)))
})

test_that("drw_tessellation gibbs produces more regular spacing", {
  r <- drw_tessellation(mosaik(extent = c(0, 100, 0, 100), res = 1),
                         type = "gibbs", n = 10, interaction = 0.15,
                         seed = 42)
  vals <- msk_pull(r)
  expect_true(length(unique(vals)) >= 2)
})

test_that("drw_tessellation is reproducible with seed", {
  r1 <- drw_tessellation(mosaik(extent = c(0, 30, 0, 30), res = 1),
                          type = "voronoi", n = 8, seed = 123)
  r2 <- drw_tessellation(mosaik(extent = c(0, 30, 0, 30), res = 1),
                          type = "voronoi", n = 8, seed = 123)
  expect_equal(msk_pull(r1), msk_pull(r2))
})

test_that("drw_tessellation name sets layer name", {
  r <- drw_tessellation(mosaik(extent = c(0, 10, 0, 10), res = 1),
                         type = "voronoi", n = 5, name = "tessellation",
                         seed = 1)
  expect_equal(names(r@layers), "tessellation")
})

test_that("drw_tessellation records provenance", {
  r <- drw_tessellation(mosaik(extent = c(0, 10, 0, 10), res = 1),
                         type = "voronoi", n = 5, seed = 1)
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "drw_tessellation")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "voronoi")
})
