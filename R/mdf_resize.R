#' Change the size of a mosaik
#'
#' Increase or decrease the size of a mosaik using nearest-neighbour
#' interpolation.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param factor [`numeric(1)`][numeric]\cr an integer larger than 1 for
#'   up-scaling and a fraction for down-scaling (see Details).
#' @return A mosaik with rescaled dimensions. Every layer is resized, since
#'   all layers share one grid, and the class labels are kept.
#' @details \code{factor} multiplies the number of cells in x and y. If
#'   \code{factor > 1}, the object is up-scaled; if \code{factor < 1}, it is
#'   down-scaled. Allowed values are integers or inverse integers (1/n).
#' @examples
#' # twice and a quarter the number of cells along each side
#' m <- mdf_resize(landscape, factor = 2)
#' msk_vis(m, .layer("cover"))
#'
#' m <- mdf_resize(landscape, factor = 1/4)
#' msk_vis(m, .layer("cover"))
#' @family operators to modify the overall object
#' @importFrom checkmate assertClass assertNumeric
#' @export

mdf_resize <- function(obj = NULL,
                       factor = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumeric(x = factor, any.missing = FALSE, len = 1)

  dims <- obj@dims

  # body ----
  if(factor > 1){
    factor <- trunc(factor)
  } else{
    denominator <- round(1/factor)
    factor <- 1/denominator
  }
  new_layers <- lapply(names(obj@layers), function(nm) {
    as.numeric(rescaleGridCpp(vals = msk_pull(obj, nm), nrow = dims[2],
                              ncol = dims[1], factorRow = factor,
                              factorCol = factor)$vals)
  })
  names(new_layers) <- names(obj@layers)
  outDims <- as.integer(round(dims * factor, 0))

  # build output ----
  .warn_measured(obj, "mdf_resize")
  .update_mosaik(obj,
                 dims = outDims,
                 layers = new_layers,
                 step = step)
}
