# `mspa` is the Soille & Vogt (2009) Fig. 1 pattern. The classification
# published alongside it lives in the internal object `.mspaPublished`
# (0 background, 1 core, 2 loop, 3 perforation, 4 islet, 5 bridge, 6 edge,
# 7 branch) rather than as a second layer, since it is a test fixture and not
# something a user of the package needs.
#
# These pin the individual steps that vignette("mspa") composes into the full
# classification, which reproduces the published result exactly (1178/1178).

test_that("the mspa dataset is present and internally consistent", {
  expect_s4_class(mspa, "mosaik")
  expect_equal(msk_names(mspa), "pattern")

  X <- msk_pull(mspa, "pattern")
  C <- .mspaPublished

  # a cell carries a class exactly when it is foreground
  expect_equal(X > 0, C > 0)
  expect_true(all(C %in% 0:7))
  expect_equal(sum(X), 710)
})

test_that("core is the distance transform thresholded at size * sqrt(2)", {
  # the reference uses a strict threshold on the euclidean distance transform,
  # with edu fixed at sqrt(2); this reproduces every core cell exactly
  m <- mdf_distance(mspa, source = "background", layer = "pattern", add = "edt")
  m <- mdf_filter(m, edt > 1 * sqrt(2),
                  add = "core")

  expect_equal(msk_pull(m, "core"),
               as.numeric(.mspaPublished == 1))
})

test_that("islets are the foreground components holding no core", {
  m <- mdf_distance(mspa, source = "background", layer = "pattern", add = "edt")
  m <- mdf_filter(m, edt > 1 * sqrt(2),
                  add = "core")
  m <- mdf_componentise(m, connectivity = 8L, layer = "pattern", add = "cc")
  m <- mdf_filter(m, pattern != 0, value = TRUE, layer = "cc") |>
      mdf_replace(old = NA, new = 0, layer = "cc")
  m <- mdf_summarise(m, by = "cc", fun = "any", layer = "core",
                 background = 0, add = "hasCore")
  m <- mdf_filter(m, pattern == 1 & hasCore == 0,
                  add = "islet")

  expect_equal(msk_pull(m, "islet"),
               as.numeric(.mspaPublished == 4))
})

test_that("bridges are connectors of at least two cells touching two cores", {
  # the size condition matters: a single cell wedged diagonally between two
  # cores touches both without running between them, and is a loop
  C <- .mspaPublished

  m <- mspa
  m@layers$published <- C          # inject the fixture for this test only
  m <- mdf_distance(m, source = "background", layer = "pattern", add = "edt")
  m <- mdf_filter(m, edt > 1 * sqrt(2),
                  add = "core")
  m <- mdf_filter(m, published == 2 | published == 5, add = "conn")
  m <- mdf_componentise(m, connectivity = 8L, layer = "conn", add = "strand")
  m <- mdf_filter(m, conn != 0, value = TRUE, layer = "strand") |>
      mdf_replace(old = NA, new = 0, layer = "strand")
  m <- mdf_componentise(m, connectivity = 8L, layer = "core", add = "coreID")
  m <- mdf_filter(m, core != 0, value = TRUE, layer = "coreID") |>
      mdf_replace(old = NA, new = 0, layer = "coreID")
  m <- mdf_summarise(m, by = "strand", fun = "n_distinct", neighbours = TRUE,
                 connectivity = 8L, layer = "coreID", background = 0,
                 add = "nCores")
  m <- mdf_summarise(m, by = "strand", fun = "n", layer = "strand",
                 background = 0, add = "strandSize")
  m <- mdf_filter(m, conn == 1 & nCores >= 2 & strandSize >= 2, add = "bridge")

  expect_equal(msk_pull(m, "bridge"), as.numeric(C == 5))
})
