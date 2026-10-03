test_that("every step records itself in the PROV shape", {
  m <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest")
  last <- m@provenance[[length(m@provenance)]]
  expect_equal(names(last), "mdf_binarise")
  e <- last[[1]]
  expect_named(e, c("wasGeneratedBy", "wasDerivedFrom", "generated",
                    "wasAssociatedWith", "atTime", "hash"))
  expect_equal(e$wasDerivedFrom, "cover")
  expect_equal(e$generated, "forest")
  expect_match(e$wasAssociatedWith, "^mosaik ")
  # defaults are recorded too, not only what was typed
  expect_true("thresh" %in% names(e$wasGeneratedBy$withArguments))
})

test_that("the input layer is resolved when none is given", {
  m <- mdf_scale(landscape, range = c(0, 1))
  e <- m@provenance[[length(m@provenance)]][[1]]
  expect_equal(e$wasDerivedFrom, "cover")
  expect_equal(e$generated, "cover")
})

test_that("the constructor records a PROV entry", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(1, 25)))
  expect_length(m@provenance, 1)
  expect_equal(names(m@provenance[[1]]), "mosaik")
})

test_that("a call inside another mosaik function is not recorded twice", {
  n <- length(landscape@provenance)
  m <- mdf_interpolate(landscape, layer = "canopy")
  expect_length(m@provenance, n + 1)
  expect_equal(names(m@provenance[[n + 1]]), "mdf_interpolate")
})

test_that("a loop is one entry, not one per iteration", {
  m <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest")
  n <- length(m@provenance)
  grow <- mdf_dilate(layer = "forest")
  m <- mdf_loop(m, grow, times = 3, layer = "forest")
  expect_length(m@provenance, n + 1)
  expect_equal(m@provenance[[n + 1]][[1]]$wasGeneratedBy$iterations, 3)
})

test_that("the history of a result replays as a recipe", {
  r <- landscape |>
    mdf_binarise(match = 47, layer = "cover", add = "forest") |>
    mdf_erode(layer = "forest", add = "core")
  again <- mdf(landscape, r)
  expect_equal(msk_pull(again, "core"), msk_pull(r, "core"))
})

test_that("msk_add takes vectors of values and records large ones by digest", {
  v <- runif(msk_ncells(landscape))
  m <- msk_add(landscape, noise = v)
  expect_equal(msk_pull(m, "noise"), v)
  e <- m@provenance[[length(m@provenance)]][[1]]
  expect_s3_class(e$wasGeneratedBy$withArguments$noise, "msk_used")
  expect_equal(e$used$noise$digest, digest::digest(v))
  # such a step cannot be repeated
  expect_error(mdf(landscape, m), "digest")
})

test_that("an argument is evaluated once, not again for the record", {
  set.seed(1)
  m <- msk_add(landscape, noise = runif(msk_ncells(landscape)))
  e <- m@provenance[[length(m@provenance)]][[1]]
  expect_equal(e$used$noise$digest, digest::digest(msk_pull(m, "noise")))
})

test_that("msk_add refuses names without a mosaik to take them from", {
  expect_error(msk_add(landscape, from = NULL, cover), "without 'from'")
  expect_error(msk_add(landscape), "nothing to add")
  expect_error(msk_add(landscape, cover = 1:3), "already present|values but")
})
