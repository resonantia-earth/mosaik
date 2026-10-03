#' Change the scale of values in a mosaik
#'
#' Change minimum and maximum value and scale all values to that range.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param range [`numeric(2)`][numeric]\cr vector of minimum and maximum value
#'   to which the values shall be scaled.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, where cell values
#'   have been rescaled.
#' @examples
#' # the canopy height scaled to 0 to 1, and to -1 to 1
#' m <- landscape |>
#'   mdf_scale(range = c(0, 1), layer = "canopy", add = "zero_one") |>
#'   mdf_scale(range = c(-1, 1), layer = "canopy", add = "minus_one_one")
#' msk_vis(m, .layer("canopy"), .layer("zero_one"), .layer("minus_one_one"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertIntegerish assertCharacter
#' @export

mdf_scale <- function(obj = NULL,
                      range = NULL,
                      layer = NULL,
                      add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

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
  .update_mosaik(obj, values = temp, step = step)
}
