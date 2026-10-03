#' Select cells based on a value range
#'
#' Transform a mosaik by setting all cells outside a lower and upper threshold
#' to \code{background}.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param lower [`numeric(1)`][numeric]\cr minimum value above which values
#'   will be selected.
#' @param upper [`numeric(1)`][numeric]\cr maximum value below which values
#'   will be selected.
#' @param background [`integerish(1)`][integer]\cr the value any cell outside the
#'   range should have.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, where cells within
#'   the range retain their value and others are set to \code{background}.
#' @examples
#' # the canopy height from 5 m, up to 20 m, and between the two
#' m <- landscape |>
#'   mdf_range(lower = 5, layer = "canopy", add = "from_5") |>
#'   mdf_range(upper = 20, layer = "canopy", add = "to_20") |>
#'   mdf_range(lower = 5, upper = 20, layer = "canopy", add = "5_to_20")
#' msk_vis(m, .layer("canopy"), .layer("from_5"), .layer("to_20"),
#'         .layer("5_to_20"))
#' @family operators to select a subset of cells
#' @importFrom checkmate assertClass assertNumeric assertIntegerish
#'   assertCharacter
#' @export

mdf_range <- function(obj = NULL,
                      lower = NULL,
                      upper = NULL,
                      background = NA,
                      layer = NULL,
                      add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumeric(lower, len = 1, any.missing = FALSE, null.ok = TRUE)
  assertNumeric(upper, len = 1, any.missing = FALSE, null.ok = TRUE)
  assertIntegerish(x = background)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  uVals <- unique(vals)

  # body ----
  temp <- vals
  if(!is.null(lower)){
    if(!min(uVals) < lower | !lower < max(uVals)){
      stop("please provide values for 'lower' within the range of the values of 'obj'.")
    }
    temp <- morphCpp(vals = temp,
                     valRows = dims[2],
                     valCols = dims[1],
                     kernel = matrix(lower, 1, 1),
                     value = uVals,
                     blend = 4,
                     merge = 12,
                     rotateKernel = FALSE,
                     strictKernel = FALSE)
  } else {
    lower <- min(uVals)
  }
  if(!is.null(upper)){
    if(!min(uVals) < upper | !upper < max(uVals)){
      stop("please provide values for 'upper' within the range of the values of 'obj'.")
    }
    temp <- morphCpp(vals = temp,
                     valRows = dims[2],
                     valCols = dims[1],
                     kernel = matrix(upper, 1, 1),
                     value = uVals,
                     blend = 3,
                     merge = 12,
                     rotateKernel = FALSE,
                     strictKernel = FALSE)
  } else {
    upper <- max(uVals)
  }
  temp[is.na(temp)] <- background

  # build output ----
  .update_mosaik(obj, values = temp, step = step)
}
