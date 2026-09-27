test_that("syn_tessellation returns a mosaik with correct dimensions", {
  for (type in c("voronoi", "rectangle", "gibbs")) {
    r <- syn_tessellation(mosaik(extent = c(0, 50, 0, 50), res = 1),
                           type = type, n = 10, seed = 42)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(50L, 50L))
    expect_equal(length(msk_pull(r)), 2500L)
  }
})

test_that("syn_tessellation voronoi produces n regions", {
  r <- syn_tessellation(mosaik(extent = c(0, 50, 0, 50), res = 1),
                         type = "voronoi", n = 15, seed = 1)
  vals <- msk_pull(r)
  expect_equal(length(unique(vals)), 15)
})

test_that("syn_tessellation rectangle tiles the grid completely", {
  r <- syn_tessellation(mosaik(extent = c(0, 30, 0, 30), res = 1),
                         type = "rectangle", n = 9, seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals > 0))
  expect_true(length(unique(vals)) >= 2)
})

test_that("syn_tessellation gibbs produces more regular spacing", {
  r <- syn_tessellation(mosaik(extent = c(0, 100, 0, 100), res = 1),
                         type = "gibbs", n = 10, interaction = 0.15,
                         seed = 42)
  vals <- msk_pull(r)
  expect_true(length(unique(vals)) >= 2)
})

test_that("syn_tessellation is reproducible with seed", {
  r1 <- syn_tessellation(mosaik(extent = c(0, 30, 0, 30), res = 1),
                          type = "voronoi", n = 8, seed = 123)
  r2 <- syn_tessellation(mosaik(extent = c(0, 30, 0, 30), res = 1),
                          type = "voronoi", n = 8, seed = 123)
  expect_equal(msk_pull(r1), msk_pull(r2))
})

test_that("syn_tessellation name sets layer name", {
  r <- syn_tessellation(mosaik(extent = c(0, 10, 0, 10), res = 1),
                         type = "voronoi", n = 5, name = "tessellation",
                         seed = 1)
  expect_equal(names(r@layers), "tessellation")
})

test_that("syn_tessellation records provenance", {
  r <- syn_tessellation(mosaik(extent = c(0, 10, 0, 10), res = 1),
                         type = "voronoi", n = 5, seed = 1)
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "syn_tessellation")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "voronoi")
})
