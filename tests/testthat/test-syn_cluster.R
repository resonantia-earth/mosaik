test_that("syn_cluster returns a mosaik with correct dimensions", {
  for (type in c("percolation", "randomCluster")) {
    r <- syn_cluster(mosaik(extent = c(0, 30, 0, 30), res = 1),
                      type = type, seed = 42)
    expect_s4_class(r, "mosaik")
    expect_equal(r@dims, c(30L, 30L))
    expect_equal(length(msk_pull(r)), 900L)
  }
})

test_that("syn_cluster percolation produces binary output", {
  r <- syn_cluster(mosaik(extent = c(0, 50, 0, 50), res = 1),
                    type = "percolation", p = 0.5, seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals %in% c(0, 1)))
})

test_that("syn_cluster percolation p controls density", {
  r_low <- syn_cluster(mosaik(extent = c(0, 100, 0, 100), res = 1),
                        type = "percolation", p = 0.2, seed = 1)
  r_high <- syn_cluster(mosaik(extent = c(0, 100, 0, 100), res = 1),
                         type = "percolation", p = 0.8, seed = 1)
  prop_low <- mean(msk_pull(r_low))
  prop_high <- mean(msk_pull(r_high))
  expect_lt(prop_low, prop_high)
  expect_equal(prop_low, 0.2, tolerance = 0.05)
  expect_equal(prop_high, 0.8, tolerance = 0.05)
})

test_that("syn_cluster randomCluster has n classes", {
  r <- syn_cluster(mosaik(extent = c(0, 50, 0, 50), res = 1),
                    type = "randomCluster", p = 0.6, n = 4, seed = 1)
  vals <- msk_pull(r)
  expect_true(all(vals > 0))
  expect_lte(length(unique(vals)), 4)
})

test_that("syn_cluster randomCluster fills all cells", {
  r <- syn_cluster(mosaik(extent = c(0, 30, 0, 30), res = 1),
                    type = "randomCluster", p = 0.5, n = 3, seed = 42)
  vals <- msk_pull(r)
  expect_true(all(vals > 0))
})

test_that("syn_cluster is reproducible with seed", {
  r1 <- syn_cluster(mosaik(extent = c(0, 30, 0, 30), res = 1),
                     type = "percolation", seed = 42)
  r2 <- syn_cluster(mosaik(extent = c(0, 30, 0, 30), res = 1),
                     type = "percolation", seed = 42)
  expect_equal(msk_pull(r1), msk_pull(r2))
})

test_that("syn_cluster name sets layer name", {
  r <- syn_cluster(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    type = "percolation", name = "cluster", seed = 1)
  expect_equal(names(r@layers), "cluster")
})

test_that("syn_cluster records provenance", {
  r <- syn_cluster(mosaik(extent = c(0, 10, 0, 10), res = 1),
                    type = "randomCluster", seed = 1)
  last_prov <- r@provenance[[length(r@provenance)]]
  expect_equal(names(last_prov)[1], "syn_cluster")
  expect_equal(last_prov[[1]]$wasGeneratedBy$withArguments$type, "randomCluster")
})
