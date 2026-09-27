#' Change the size of a mosaik
#'
#' Increase or decrease the size of a mosaik using nearest-neighbour
#' interpolation.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param factor [numeric(1)][numeric]\cr an integer larger than 1 for
#'   up-scaling and a fraction for down-scaling (see Details).
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @return A mosaik with rescaled dimensions.
#' @details \code{factor} multiplies the number of cells in x and y. If
#'   \code{factor > 1}, the object is up-scaled; if \code{factor < 1}, it is
#'   down-scaled. Allowed values are integers or inverse integers (1/n).
#' @examples
#' mdf_resize(landscape, factor = 2)
#' mdf_resize(landscape, factor = 1/4)
#' @family operators to modify the overall object
#' @importFrom checkmate assertClass assertNumeric assertCharacter
#' @export

mdf_resize <- function(obj = NULL,
                       factor = NULL,
                       layer = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumeric(x = factor, any.missing = FALSE, len = 1)
  assertCharacter(x = layer, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims

  # body ----
  if(factor > 1){
    factor <- trunc(factor)
  } else{
    denominator <- round(1/factor)
    factor <- 1/denominator
  }
  result <- rescaleGridCpp(vals = vals,
                           nrow = dims[2],
                           ncol = dims[1],
                           factorRow = factor,
                           factorCol = factor)
  outDims <- as.integer(round(dims * factor, 0))

  # build output ----
  prov <- msk_prov("mdf_resize", list(factor = factor, layer = layer))
  new_mosaik(
    extent     = obj@extent,
    dims       = outDims,
    layers     = stats::setNames(list(as.numeric(result$vals)), layer),
    categories = list(),
    patches    = list(),
    global     = list(),
    crs        = obj@crs,
    provenance = c(obj@provenance, list(prov))
  )
}
