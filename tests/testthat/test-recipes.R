# Every recipe of the library runs, and its check holds: each recipe ends its
# check with a hidden stopifnot(), which fails the test when the result no
# longer is what the recipe's text predicts.

recipes <- list.files(system.file("recipes", package = "mosaik"),
                      pattern = "[.]Rmd$", full.names = TRUE)
vocabulary <- yaml::read_yaml(system.file("recipes", "facets.yml",
                                          package = "mosaik"))
objects <- vocabulary$objects
facets <- vocabulary[setdiff(names(vocabulary), c("objects", "kinds"))]

test_that("the recipe library is there", {
  expect_gt(length(recipes), 0)
  # an output names the objects it is derived from, and those exist
  from <- unlist(lapply(objects, `[[`, "from"))
  expect_true(all(from %in% names(objects)))
  side <- vapply(objects, `[[`, "", "side")
  expect_true(all(side %in% c("input", "state", "output")))
  # an input map is of a known kind
  kind <- unlist(lapply(objects[side == "input"], `[[`, "kind"))
  expect_true(all(kind %in% names(vocabulary$kinds)))
  expect_length(kind, sum(side == "input"))
})

for (rmd in recipes) {
  test_that(paste("recipe", basename(rmd), "names known facets"), {
    head <- rmarkdown::yaml_front_matter(rmd)
    for (f in names(facets)) {
      expect_true(head[[f]] %in% names(facets[[f]]$values),
                  info = paste(f, "=", head[[f]]))
    }
    expect_true(length(head$uses) > 0)
    expect_true(length(head$result) > 0)
    expect_true(all(head$uses %in% names(objects)),
                info = paste(setdiff(head$uses, names(objects)), collapse = ", "))
    expect_true(nzchar(head$title))
    expect_true(nzchar(head$summary))
  })

  test_that(paste("recipe", basename(rmd), "runs and its check holds"), {
    skip_on_cran()
    code <- tempfile(fileext = ".R")
    knitr::purl(rmd, output = code, quiet = TRUE, documentation = 0)
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    env <- new.env(parent = globalenv())
    expect_no_error(sys.source(code, envir = env))
    # every result the overview reads exists after the recipe has run
    for (r in rmarkdown::yaml_front_matter(rmd)$result) {
      expect_true(exists(r$object, envir = env, inherits = FALSE), info = r$object)
      expect_true(r$layer %in% msk_names(get(r$object, envir = env)),
                  info = paste(r$object, r$layer))
    }
  })
}
