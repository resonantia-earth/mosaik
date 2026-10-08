#' Perimeter length
#'
#' Calculate the length of the boundary of each class of a layer and store it
#' in the class table of that layer.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param unit [`character(1)`][character]\cr \code{"map"} (default, in map
#'   units) or \code{"cells"} (number of cell edges).
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with \code{perimeter} added to the class table of
#'   \code{layer} (see \code{\link{msk_table}}).
#' @details The boundary of a class runs between its cells and every
#'   neighbouring cell that is not of the class, including cells that are
#'   \code{NA}. On a layer of patch numbers from \code{\link{mdf_componentise}},
#'   the cells outside every patch are \code{NA}, so each patch has its full
#'   outline. Edges along the map border are not counted: the map ends there,
#'   the class may not. The perimeter of the whole layer is
#'   \code{sum(perimeter.all)} in \code{\link{msr}}.
#' @examples
#' # the perimeter of each land cover class
#' m <- msr_perimeter(landscape, layer = "cover")
#' msk_table(m, layer = "cover")$perimeter
#'
#' # the perimeter of each forest patch
#' m <- mdf_filter(landscape, cover == 47, add = "forest") |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
#'   msr_perimeter(layer = "patch")
#' msk_table(m, layer = "patch")$perimeter
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_perimeter <- function(obj = NULL, unit = "map", layer = NULL){

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

  # an NA cell is not of the class, so the edge to it counts: NA becomes a
  # value below every class, whose own edges are dropped again
  outside <- if (all(is.na(vals))) 0 else min(vals, na.rm = TRUE) - 1
  vals[is.na(vals)] <- outside

  values <- countCellEdgesCpp(vals = as.numeric(vals), nrow = dims[2],
                              ncol = dims[1])
  values <- values[values$value != outside, ]

  if(unit == "map"){
    edges <- values$edgesX * theRes[2] + values$edgesY * theRes[1]
  } else {
    edges <- values$edgesX + values$edgesY
  }

  obj <- .store_class(obj, layer, "perimeter", as.integer(values$value), edges)

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
