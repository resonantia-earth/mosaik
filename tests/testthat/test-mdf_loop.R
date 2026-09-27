seeded <- function(){
  # two patches, only the left one holds a seed cell
  g <- matrix(0L, 7, 11)
  g[2:6, 2:5]  <- 1L
  g[2:6, 7:10] <- 1L
  s <- matrix(0L, 7, 11); s[4, 3] <- 1L
  mosaik(extent = c(0, 11, 0, 7), res = 1,
         vals = list(mask = as.numeric(t(g)), seed = as.numeric(t(s))))
}
growRecipe <- function(){
  mdf_dilate(struct = msk_struct("square", width = 3), layer = "seed") |>
    mdf_filter(seed == 1 & mask == 1, value = TRUE, background = 0,
               add = "seed")
}

test_that("stable runs to a fixpoint", {
  # the seed fills its own patch entirely and never reaches the other one
  o <- mdf_loop(seeded(), growRecipe(), stable = TRUE, layer = "seed")
  s <- matrix(msk_pull(o, "seed"), nrow = 7, byrow = TRUE)

  expect_equal(sum(s), 20)              # the left patch, 5 x 4
  expect_true(all(s[2:6, 2:5] == 1))
  expect_true(all(s[, 7:10] == 0))      # the unseeded patch stays empty
})

test_that("times runs a fixed number of iterations", {
  one <- mdf_loop(seeded(), growRecipe(), times = 1)
  two <- mdf_loop(seeded(), growRecipe(), times = 2)

  expect_equal(sum(msk_pull(one, "seed")), 9)    # 3x3 around the seed
  expect_lt(sum(msk_pull(one, "seed")), sum(msk_pull(two, "seed")))
})

test_that("until stops on a condition over the layers", {
  peel <- mdf_erode(struct = msk_struct("square", width = 3), layer = "mask")
  o <- mdf_loop(seeded(), peel, until = sum(mask, na.rm = TRUE) == 0)
  expect_equal(sum(msk_pull(o, "mask"), na.rm = TRUE), 0)
})

test_that("times caps a condition that would not settle", {
  grow <- mdf_dilate(struct = msk_struct("square", width = 3), layer = "seed")
  # unmasked dilation never stabilises within the grid, so times must stop it
  o <- mdf_loop(seeded(), grow, times = 2, stable = TRUE, layer = "seed")
  expect_gt(sum(msk_pull(o, "seed")), 9)
})

test_that("an endless loop is refused", {
  expect_error(mdf_loop(seeded(), growRecipe()), "stopping condition")
})

test_that("mdf_loop is itself recordable and replays", {
  r <- mdf_loop(obj = NULL, growRecipe(), stable = TRUE, layer = "seed")
  o <- mdf(seeded(), r)
  expect_equal(sum(msk_pull(o, "seed")), 20)
})

test_that("the number of iterations is recorded in the provenance", {
  o <- mdf_loop(seeded(), growRecipe(), stable = TRUE, layer = "seed")
  last <- o@provenance[[length(o@provenance)]]
  expect_equal(names(o@provenance)[length(o@provenance)], "mdf_loop")
  expect_gt(last$wasGeneratedBy$withArguments$iterations, 1)
})
