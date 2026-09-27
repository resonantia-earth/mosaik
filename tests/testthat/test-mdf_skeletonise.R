test_that("zhangSuen thins a bar to a one-cell line", {
  g <- matrix(0L, nrow = 9, ncol = 9)
  g[4:6, 2:8] <- 1L
  m <- mosaik(extent = c(0, 9, 0, 9), res = 1,
              vals = list(X = as.numeric(t(g))))
  s <- msk_pull(mdf_skeletonise(m, background = 0, layer = "X", add = "s"), "s")

  expect_true(all(s %in% c(0, 1)))
  expect_lt(sum(s), sum(g))          # it thinned
  expect_gt(sum(s), 0)               # but kept the bar
})

test_that("homotopic thinning collapses a simply-connected blob to a point", {
  # a solid block is topologically a point, and homotopic thinning preserves
  # topology -- so unanchored it must reduce to a single cell. This is the
  # defining difference from zhangSuen, which protects line tips.
  g <- matrix(0L, nrow = 9, ncol = 9)
  g[2:8, 2:8] <- 1L
  m <- mosaik(extent = c(0, 9, 0, 9), res = 1,
              vals = list(X = as.numeric(t(g))))
  s <- msk_pull(mdf_skeletonise(m, background = 0, layer = "X", add = "s",
                                method = "homotopic"), "s")

  expect_equal(sum(s), 1)
})

test_that("homotopic thinning keeps a ring, which is not simply connected", {
  g <- matrix(0L, nrow = 9, ncol = 9)
  g[2:8, 2:8] <- 1L
  g[4:6, 4:6] <- 0L                  # punch a hole -> an annulus
  m <- mosaik(extent = c(0, 9, 0, 9), res = 1,
              vals = list(X = as.numeric(t(g))))
  s <- msk_pull(mdf_skeletonise(m, background = 0, layer = "X", add = "s",
                                method = "homotopic"), "s")

  # the loop around the hole cannot be removed without changing topology
  expect_gt(sum(s), 1)
})

test_that("homotopic thinning is independent of the scan order", {
  # The defining property: the result must be invariant under transformations
  # that change which cell a scan reaches first. Deleting simple points one at
  # a time does NOT give this -- two adjacent simple pixels can each be
  # removable alone but not together, so whichever is visited first wins.
  set.seed(42)
  g <- matrix(0L, nrow = 16, ncol = 16)
  g[3:13, 3:13] <- rbinom(11 * 11, 1, 0.75)
  m <- mosaik(extent = c(0, 16, 0, 16), res = 1,
              vals = list(X = as.numeric(t(g))))
  s <- matrix(msk_pull(mdf_skeletonise(m, background = 0, layer = "X",
                                       add = "s", method = "homotopic"), "s"),
              nrow = 16, byrow = TRUE)

  # flipping the grid must flip the skeleton, not change it
  for(flip in list(function(x) x[nrow(x):1, ],
                   function(x) x[, ncol(x):1],
                   function(x) t(x))){
    gf <- flip(g)
    mf <- mosaik(extent = c(0, 16, 0, 16), res = 1,
                 vals = list(X = as.numeric(t(gf))))
    sf <- matrix(msk_pull(mdf_skeletonise(mf, background = 0, layer = "X",
                                          add = "s", method = "homotopic"), "s"),
                 nrow = 16, byrow = TRUE)
    expect_equal(sf, flip(s))
  }
})

test_that("anchored homotopic thinning retracts dead ends but keeps connectors", {
  # two blocks joined by a neck, plus a dead-end arm hanging off the left block.
  # Anchored on the blocks: the neck is needed to keep them connected as one
  # component and survives; the arm is not needed and retracts entirely. That
  # retraction is what MSPA relies on to tell a connector from a branch.
  g <- matrix(0L, nrow = 9, ncol = 17)
  g[3:7, 2:6]   <- 1L                # left block
  g[3:7, 12:16] <- 1L                # right block
  g[5, 7:11]    <- 1L                # neck joining them
  g[2, 3:5]     <- 1L                # dead-end arm on top of the left block
  a <- matrix(0L, nrow = 9, ncol = 17)
  a[4:6, 3:5]   <- 1L
  a[4:6, 13:15] <- 1L

  m <- mosaik(extent = c(0, 17, 0, 9), res = 1,
              vals = list(X = as.numeric(t(g)), anchor = as.numeric(t(a))))

  s <- matrix(msk_pull(mdf_skeletonise(m, background = 0, anchor = "anchor",
                                       layer = "X", add = "s",
                                       method = "homotopic"), "s"),
              nrow = 9, byrow = TRUE)

  expect_true(all(s[5, 7:11] == 1))  # the connector survives
  expect_true(all(s[2, 3:5] == 0))   # the dead-end arm is gone
})
