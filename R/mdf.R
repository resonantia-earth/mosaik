#' Apply a recipe to a mosaik
#'
#' Run a recorded sequence of \code{mdf_*} and \code{msr_*} calls on a real
#' mosaik. A recipe is built by calling those functions with \code{obj = NULL}
#' (or chaining them from an empty mosaik): each call records itself instead of
#' executing, and the accumulated \code{@provenance} is the recipe. \code{mdf()}
#' replays that sequence on \code{obj}, threading the result of each step into
#' the next, so an algorithm can be defined once and applied to many rasters.
#'
#' Because every call records itself the same way, the history of any mosaik
#' is a recipe too: \code{mdf(new, result)} repeats on \code{new} what was done
#' to \code{result}. The entry that created \code{result} is skipped, and a
#' step whose input was too large to record (a whole vector of values, another
#' mosaik) cannot be repeated and stops with an error.
#'
#' A recipe repeats operations, not data. Every argument is recorded as it was
#' given, so a value read off one map, such as a list of patch numbers passed
#' to \code{\link{mdf_replace}}, is a constant in the recipe and stays the same
#' on every other map. A recipe built with \code{obj = NULL} has no data to read
#' values from, so it is the way to be sure a recipe runs on other inputs as it
#' is.
#'
#' Steps that read an earlier intermediate (e.g. the \code{by} of
#' \code{\link{mdf_summarise}}) reference it by the layer name an earlier step wrote
#' via \code{add}. The chain therefore stays linear in execution while still
#' expressing branch-and-recombine algorithms through named layers.
#'
#' @param obj [`mosaik`]\cr the mosaik to apply the recipe to.
#' @param recipe [`mosaik`]\cr a recipe shell, or any mosaik, whose
#'   \code{@provenance} holds the steps to repeat.
#' @return The mosaik after all recipe steps have been applied.
#' @examples
#' # define once, without data: the core of class 47, one cell in from its edge
#' core <- mdf_filter(expr = cover == 47, add = "forest") |>
#'   mdf_erode(layer = "forest", add = "core")
#'
#' # apply to any raster with a 'cover' layer
#' result <- mdf(landscape, core)
#' msk_vis(result, .layer("cover"), .layer("forest"), .layer("core"))
#'
#' rotated <- mdf(mdf_rotate(landscape, angle = 90), core)
#' msk_vis(rotated, .layer("cover"), .layer("forest"), .layer("core"))
#' @importFrom checkmate assertClass
#' @export

mdf <- function(obj, recipe){

  assertClass(x = obj, classes = "mosaik")
  assertClass(x = recipe, classes = "mosaik")

  # every recorded call except the one that created the recipe's object
  # each history element is list(<function> = <entry>)
  prov  <- recipe@provenance
  prov  <- prov[vapply(prov, function(e) is.list(e) && !is.null(e[[1]]$wasGeneratedBy),
                       logical(1))]
  steps <- lapply(prov, `[[`, 1)
  names(steps) <- vapply(prov, names, character(1))
  steps <- steps[names(steps) != "mosaik"]
  if (length(steps) == 0) {
    stop("'recipe' contains no recorded steps. Build it by calling mdf_* or msr_* ",
         "with obj = NULL.", call. = FALSE)
  }
  lost <- vapply(steps, function(e) any(vapply(.prov_args(e), inherits,
                                               logical(1), "msk_used")),
                 logical(1))
  if (any(lost)) {
    stop("step(s) ", paste(which(lost), collapse = ", "), " (",
         paste(unique(names(steps)[lost]), collapse = ", "), ") recorded an ",
         "input only by its digest and cannot be repeated.", call. = FALSE)
  }

  for (i in seq_along(steps)) {
    fn   <- names(steps)[i]
    args <- .prov_args(steps[[i]])
    obj  <- do.call(fn, c(list(obj = obj), args))
  }

  obj
}
