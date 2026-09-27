#' Change the scale of values in a mosaik
#'
#' Change minimum and maximum value and scale all values to that range.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param range [numeric(2)][numeric]\cr vector of minimum and maximum value
#'   to which the values shall be scaled.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, where cell values
#'   have been rescaled.
#' @examples
#' mdf_scale(landscape, range = c(0, 1), layer = "intensity")
#' mdf_scale(landscape, range = c(-1, 1), layer = "intensity")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertIntegerish assertCharacter
#' @export

mdf_scale <- function(obj = NULL,
                      range = NULL,
                      layer = NULL,
                      add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertIntegerish(x = range, len = 2)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  minVal <- min(vals, na.rm = TRUE)
  maxVal <- max(vals, na.rm = TRUE)
  temp <- (vals - minVal) * (range[2] - range[1]) / (maxVal - minVal) + range[1]

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_scale", list(range = range, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
