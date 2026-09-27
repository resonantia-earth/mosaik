test_that(".themeDefaults holds the default settings", {
  th <- mosaik:::.themeDefaults
  expect_type(th, "list")
  expect_true(th$panel.title.plot)
  expect_true(th$figure.title.plot)
  expect_true(th$box.plot)
  expect_false(any(duplicated(names(th))))
})

test_that(".complete_theme fills in from the defaults", {
  th <- mosaik:::.complete_theme(.theme(panel.title.plot = FALSE))
  expect_false(th$panel.title.plot)
  expect_equal(th$panel.title.fontsize,
               mosaik:::.themeDefaults$panel.title.fontsize)
  expect_equal(length(th), length(mosaik:::.themeDefaults))
  # a bare list is not a theme
  expect_error(mosaik:::.complete_theme(list(panel.title.plot = FALSE)),
               "must be built with .theme")
})

test_that(".theme validates and extends", {
  th <- .theme(legend.plot = FALSE)
  expect_s3_class(th, "mskTheme")
  expect_false(th$legend.plot)

  # extending keeps the earlier settings
  ext <- .theme(th, panel.title.fontsize = 14)
  expect_false(ext$legend.plot)
  expect_equal(ext$panel.title.fontsize, 14)
  expect_error(.theme(list(a = 1)), "must be a theme built with")

  expect_error(.theme(nonsense = 1), "unknown theme")
  # a near miss gets a suggestion
  expect_error(.theme(panel.title.fontsze = 1), "did you mean")
  expect_error(.theme(TRUE), "must be a theme built with")
})

test_that(".makeScale bins continuous layers and keeps categories", {
  expect_length(mosaik:::.makeScale(c(1, 2, 3, 4, 5)), 5)
  expect_length(mosaik:::.makeScale(runif(5000), bins = 256), 256)
  expect_length(mosaik:::.makeScale(runif(5000), bins = 10), 10)
  expect_equal(range(mosaik:::.makeScale(runif(500), limits = c(0, 10))),
               c(0, 10))
  expect_equal(mosaik:::.makeScale(rep(7, 10)), 7)
  expect_length(mosaik:::.makeScale(numeric(0)), 0)
})

test_that(".resolve_colours handles ramps, palettes and category maps", {
  expect_length(mosaik:::.resolve_colours(c("black", "white"), 5), 5)
  expect_length(mosaik:::.resolve_colours("viridis", 7), 7)
  # a named vector is a category map and comes back untouched
  cm <- c(a = "red", b = "blue")
  expect_equal(mosaik:::.resolve_colours(cm, 2), cm)
  expect_error(mosaik:::.resolve_colours("notapalette", 5), "neither a colour")
})

test_that(".hillshade is directional and bounded", {
  m <- mosaik(extent = c(0, 1000, 0, 1000), res = 10)
  n <- prod(msk_dims(m))
  # a flat surface must return exactly cos(zenith)
  m@layers$flat <- rep(100, n)
  expect_equal(unique(mosaik:::.hillshade(m, "flat")), cos(45 * pi / 180))

  # a cone has every aspect, so NW-facing flanks must be lit and SE-facing
  # flanks shaded under a north-west light
  nr <- msk_dims(m)[2]; nc <- msk_dims(m)[1]
  rr <- row(matrix(0, nr, nc)); cc <- col(matrix(0, nr, nc))
  d <- sqrt((rr - nr / 2)^2 + (cc - nc / 2)^2)
  m@layers$cone <- as.numeric(t(pmax(0, 300 * (1 - d / (nr / 2)))))
  sh <- mosaik:::.hillshade(m, "cone", azimuth = 315)
  sa <- mosaik:::slopeAspectCpp(as.numeric(msk_pull(m, "cone")), nr, nc, 10, 10)
  on <- as.numeric(msk_pull(m, "cone")) > 1
  d2r <- pi / 180
  faceNW <- on & sa$aspect > 270 * d2r & sa$aspect < 360 * d2r
  faceSE <- on & sa$aspect > 90 * d2r & sa$aspect < 180 * d2r
  expect_gt(mean(sh[faceNW]), mean(sh[faceSE]))
  expect_true(all(sh >= 0 & sh <= 1))
})

test_that(".hillshade averages several azimuths", {
  m <- mosaik(extent = c(0, 1000, 0, 1000), res = 10)
  nr <- msk_dims(m)[2]; nc <- msk_dims(m)[1]
  rr <- row(matrix(0, nr, nc)); cc <- col(matrix(0, nr, nc))
  d <- sqrt((rr - nr / 2)^2 + (cc - nc / 2)^2)
  m@layers$cone <- as.numeric(t(pmax(0, 300 * (1 - d / (nr / 2)))))
  on <- as.numeric(msk_pull(m, "cone")) > 1
  one <- mosaik:::.hillshade(m, "cone", azimuth = 315)
  many <- mosaik:::.hillshade(m, "cone", azimuth = c(225, 270, 315, 360))
  # spreading the light over four directions lifts the darkest faces, so the
  # shading spans a narrower range than a single light does
  expect_gt(min(many[on]), min(one[on]))
  expect_lt(diff(range(many[on])), diff(range(one[on])))
  expect_true(all(many >= 0 & many <= 1))
})

test_that(".apply_hillshade darkens and preserves NA", {
  cols <- c("#808080", "#808080", NA)
  expect_equal(mosaik:::.apply_hillshade(cols, c(1, 1, 1), 1)[1], "#808080FF")
  expect_equal(mosaik:::.apply_hillshade(cols, c(0.5, 0.5, 0.5), 1)[1],
               "#404040FF")
  # intensity 0 is a no-op
  expect_equal(mosaik:::.apply_hillshade(cols, c(0.2, 0.2, 0.2), 0)[1],
               "#808080FF")
  expect_true(is.na(mosaik:::.apply_hillshade(cols, c(1, 1, 1), 1)[3]))
})

test_that(".layer validates its specification", {
  s <- .layer("cover", panel = "a", colours = "viridis")
  expect_s3_class(s, "mskLayer")
  expect_equal(s$layer, "cover")
  expect_error(.layer("a", hillshade = list(azimuth = 10)), "needs a 'layer'")
  expect_error(.layer("a", hillshade = list(layer = "d", azimut = 1)),
               "subset")
  expect_error(.layer("a", hillshade = list(layer = "d", intensity = 5)))
  expect_error(.layer("a", limits = c(10, 1)), "sorted")
})

test_that("msk_vis renders without error", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(cover = sample(1:5, 100, replace = TRUE)))
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  expect_silent(msk_vis(m))
  expect_silent(msk_vis(m, .layer("cover")))
  dev.off()
  expect_true(file.exists(f))
  unlink(f)
})

test_that("msk_vis draws one panel per layer by default", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(a = rep(1, 25), b = rep(2, 25)))
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  expect_silent(msk_vis(m))
  dev.off()
  unlink(f)
})

test_that("msk_vis composites layers sharing a panel", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(dem = as.numeric(1:100),
                          riv = c(rep(NA, 90), rep(1, 10))))
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  expect_silent(msk_vis(m,
                        .layer("dem", panel = "one", colours = "terrain"),
                        .layer("riv", panel = "one", colours = "blue")))
  dev.off()
  unlink(f)
})

test_that("msk_vis rejects non-specs and missing layers", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(a = rep(1, 25)))
  expect_error(msk_vis(m, "cover"), "must be a layer specification")
  expect_error(msk_vis(m, .layer("nope")), "not in this mosaik")
})

test_that("msk_vis accepts a .theme() and rejects typos", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1, vals = list(cover = 1:25))
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  expect_silent(msk_vis(m, theme = .theme(panel.title.plot = FALSE,
                                          legend.plot = FALSE)))
  expect_error(msk_vis(m, theme = .theme(titel.plot = FALSE)), "unknown theme")
  expect_error(msk_vis(m, theme = list(panel.title.plot = FALSE)),
               "must be built with .theme")
  dev.off()
  unlink(f)
})

test_that("shared_scale pools a layer's range across panels", {
  m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
              vals = list(a = as.numeric(1:100), b = as.numeric(1:100) * 3))
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  expect_silent(msk_vis(m, .layer("a"), .layer("b"), shared_scale = TRUE))
  dev.off()
  unlink(f)
})

test_that("a named colour vector maps categories explicitly", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(c(1, 2), length.out = 25)))
  m@categories$cover <- list(gid = c(1, 2), val = c("forest", "crop"))
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  expect_silent(msk_vis(m, .layer("cover",
                                  colours = c(forest = "darkgreen",
                                              crop = "khaki"))))
  # an incomplete map is an error, not a silent fallback
  expect_error(msk_vis(m, .layer("cover", colours = c(forest = "darkgreen"))),
               "no colour given")
  dev.off()
  unlink(f)
})

test_that("msk_vis with trace prints provenance", {
  m <- mosaik(extent = c(0, 5, 0, 5), res = 1,
              vals = list(cover = rep(1, 25)))
  m2 <- mdf_offset(obj = m, value = 5)
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  expect_message(msk_vis(m2, trace = TRUE), "history")
  dev.off()
  unlink(f)
})
