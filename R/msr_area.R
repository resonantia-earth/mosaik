#' Area of classes
#'
#' Calculate the area of each class of a layer and store it in the class table
#' of that layer.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param unit [`character(1)`][character]\cr \code{"cells"} (default, number of
#'   cells) or \code{"map"} (in map units).
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with \code{area} added to the class table of
#'   \code{layer} (see \code{\link{msk_categories}}).
#' @details The classes of a layer can be land cover classes or groups, such
#'   as the patches numbered by \code{\link{mdf_componentise}}: on a layer of
#'   patch numbers, each class is a patch. The area of the whole layer is the
#'   sum over its classes, \code{sum(area.all)} in \code{\link{msr}}.
#' @examples
#' # the area of each land cover class
#' m <- msr_area(landscape, layer = "cover")
#' msk_categories(m, layer = "cover")$area
#'
#' # the area of each forest patch
#' m <- mdf_filter(landscape, cover == 47, add = "forest") |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
#'   msr_area(layer = "patch")
#' msk_categories(m, layer = "patch")$area
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_area <- function(obj = NULL, unit = "cells", layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = unit, choices = c("cells", "map"))
  assertCharacter(x = layer, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  theRes <- msk_res(obj)

  values <- countCellValuesCpp(vals = vals, nrow = dims[2], ncol = dims[1])
  area <- values$cells
  if(unit == "map") area <- area * theRes[1] * theRes[2]

  obj <- .store_class(obj, layer, "area", as.integer(values$value), area)

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
