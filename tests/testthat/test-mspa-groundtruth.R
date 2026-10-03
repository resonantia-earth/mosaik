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

test_that("the composed classification reproduces the published one cell for cell", {
  # the whole of MSPA from operators: core, islet, edge and perforation by
  # peeling, connectors by homotopic thinning anchored on the core, then
  # bridge, loop and branch; written with size 1 and edu sqrt(2), as the
  # reference implementation fixes them
  size <- 1; edu <- sqrt(2)
  m <- mdf_distance(mspa, source = "background", layer = "pattern", add = "edt") |>
    mdf_filter(edt > size * edu, add = "core") |>
    mdf_componentise(connectivity = 8L, layer = "pattern", add = "cc")
  m <- mdf_filter(m, pattern != 0, value = TRUE, layer = "cc") |>
    mdf_replace(old = NA, new = 0, layer = "cc") |>
    mdf_summarise(by = "cc", fun = "any", layer = "core", background = 0,
                  add = "hasCore") |>
    mdf_filter(pattern == 1 & hasCore == 0, add = "patch")

  # edge and perforation: peel the filled core level by level
  m <- mdf_distance(m, source = "foreground", layer = "core", add = "dcore") |>
    mdf_filter(core == 1 | (dcore > 0 & dcore <= size * edu), add = "opening") |>
    mdf_offset(value = 0, layer = "core", add = "work") |>
    mdf_fill(connectivity = 4L, layer = "work", add = "corefill") |>
    mdf_filter(corefill == 1 & work == 0, add = "holes") |>
    mdf_filter(pattern < 0, add = "edges")
  # a recipe keeps its predicates as expressions and evaluates them when it is
  # replayed, where local variables are out of reach, so the width is written
  # out: size * edu = sqrt(2)
  peel <- mdf_distance(source = "foreground", layer = "corefill", add = "dfill") |>
    mdf_filter(edges == 1 | (corefill == 0 & dfill > 0 & dfill <= sqrt(2)),
               add = "edges") |>
    mdf_fill(connectivity = 8L, layer = "holes", add = "i0") |>
    mdf_filter(corefill == 1 & i0 == 0, add = "corefill") |>
    mdf_filter(work == 1 & corefill == 0, add = "work") |>
    mdf_fill(connectivity = 4L, layer = "work", add = "corefill") |>
    mdf_filter(corefill == 1 & work == 0, add = "holes")
  m <- mdf_loop(m, peel, until = sum(corefill, na.rm = TRUE) == 0) |>
    mdf_filter(edges == 1 & core == 0, add = "edges") |>
    mdf_filter(core == 0 & dcore > 0 & dcore <= size * edu & edges == 0,
               add = "perf")

  # connectors: what survives thinning anchored on the core
  m <- mdf_filter(m, pattern == 1 & core == 0 & patch == 0 & perf == 0 &
                    edges == 0, add = "residues") |>
    mdf_filter(opening == 1 | residues == 1, add = "skin") |>
    mdf_skeletonise(anchor = "core", background = 0, layer = "skin", add = "sk",
                    method = "homotopic") |>
    mdf_filter(sk == 1 & core == 0, add = "sk") |>
    mdf_filter((opening == 1 & core == 0) | residues == 1, add = "i0") |>
    mdf_offset(value = 0, layer = "sk", add = "rec")
  grow <- mdf_dilate(struct = msk_struct("square", width = 3), layer = "rec") |>
    mdf_filter(rec == 1 & i0 == 1, add = "rec")
  m <- mdf_loop(m, grow, stable = TRUE, layer = "rec") |>
    mdf_filter(sk == 1 & i0 == 1 & rec == 1, add = "connector") |>
    mdf_componentise(connectivity = 8L, layer = "connector", add = "strand")
  m <- mdf_filter(m, connector != 0, value = TRUE, layer = "strand") |>
    mdf_replace(old = NA, new = 0, layer = "strand") |>
    mdf_componentise(connectivity = 8L, layer = "core", add = "coreID")
  m <- mdf_filter(m, core != 0, value = TRUE, layer = "coreID") |>
    mdf_replace(old = NA, new = 0, layer = "coreID") |>
    mdf_summarise(by = "strand", fun = "n_distinct", neighbours = TRUE,
                  connectivity = 8L, layer = "coreID", background = 0,
                  add = "nCores") |>
    mdf_summarise(by = "strand", fun = "n", layer = "strand", background = 0,
                  add = "strandSize") |>
    mdf_filter(connector == 1 & nCores >= 2 & strandSize >= 2, add = "bridge") |>
    mdf_filter(connector == 1 & bridge == 0, add = "loop") |>
    mdf_filter(pattern == 1 & core == 0 & patch == 0 & perf == 0 & edges == 0 &
                 bridge == 0 & loop == 0, add = "branch")

  # one map of codes; a connector takes precedence over a boundary
  codes <- c(branch = 7, edges = 6, perf = 3, patch = 4, core = 1, bridge = 5,
             loop = 2)
  cls <- rep(0, length(.mspaPublished))
  for (nm in names(codes)) cls[msk_pull(m, nm) == 1] <- codes[[nm]]
  cls[msk_pull(m, "pattern") != 1] <- 0

  expect_equal(sum(cls == .mspaPublished), 1178)
})
