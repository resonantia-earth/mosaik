#' Apply a recipe to a mosaik
#'
#' Run a recorded sequence of \code{mdf_*} calls on a real mosaik. A recipe is built by calling those functions with \code{obj = NULL}
#' (or chaining them from an empty mosaik): each call records itself instead of
#' executing, and the accumulated \code{@provenance} is the recipe. \code{mdf()}
#' replays that sequence on \code{obj}, threading the result of each step into
#' the next, so an algorithm can be defined once and applied to many rasters.
#'
#' Steps that read an earlier intermediate (e.g. the \code{by} of
#' \code{\link{mdf_mask}}) reference it by the layer name an earlier step wrote
#' via \code{add}. The chain therefore stays linear in execution while still
#' expressing branch-and-recombine algorithms through named layers.
#'
#' @param obj [mosaik]\cr the mosaik to apply the recipe to.
#' @param recipe [mosaik]\cr a recipe shell whose \code{@provenance} holds the
#'   recorded steps (see Details).
#' @return The mosaik after all recipe steps have been applied.
#' @examples
#' # define once, without data: the core of class 47, one cell in from its edge
#' core <- mdf_binarise(match = 47, layer = "cover", add = "forest") |>
#'   mdf_erode(layer = "forest", add = "core")
#'
#' # apply to any raster with a 'cover' layer
#' result <- mdf(landscape, core)
#' msk_names(result)
#'
#' rotated <- mdf(mdf_rotate(landscape, angle = 90), core)
#' @importFrom checkmate assertClass
#' @export

mdf <- function(obj, recipe){

  assertClass(x = obj, classes = "mosaik")
  assertClass(x = recipe, classes = "mosaik")

  prov  <- recipe@provenance
  steps <- if (is.list(prov)) {
    Filter(function(e) is.list(e) && isTRUE(e$step), prov)
  } else {
    list()
  }
  if (length(steps) == 0) {
    stop("'recipe' contains no recorded steps. Build it by calling mdf_* ",
         "with obj = NULL.", call. = FALSE)
  }

  for (i in seq_along(steps)) {
    fn   <- names(steps)[i]
    args <- .prov_args(steps[[i]])
    obj  <- do.call(fn, c(list(obj = obj), args))
  }

  obj
}
