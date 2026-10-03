#' Morphologically dilate a mosaik
#'
#' The morphological operation 'dilate' increases the value of cells that match
#' a structuring element.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param struct [`struct(1)`][struct]\cr the structuring element; see
#'   \code{\link{msk_struct}} for details.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, where a dilation has
#'   been performed.
#' @examples
#' # the forest grown by one cell, and by two cells with a 5 x 5 square
#' m <- landscape |>
#'   mdf_binarise(match = 47, layer = "cover", add = "forest") |>
#'   mdf_dilate(layer = "forest", add = "dilated") |>
#'   mdf_dilate(struct = msk_struct("square", width = 5), layer = "forest",
#'              add = "dilated_5")
#' msk_vis(m, .layer("forest"), .layer("dilated"), .layer("dilated_5"))
#' @family operators to morphologically modify a raster
#' @importFrom checkmate assertClass assertCharacter
#' @export

mdf_dilate <- function(obj = NULL,
                       struct = NULL,
                       layer = NULL,
                       add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertClass(x = struct, classes = "struct", null.ok = TRUE)
  if(is.null(struct)) struct <- msk_struct(type = "diamond", width = 3, height = 3)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  uVals <- unique(vals)

  # body ----
  blend <- 1 # blendIdentity
  if(!isBinaryCpp(vals = vals)){
    values <- uVals[uVals != 0]
    if(!isBinaryCpp(vals = as.numeric(struct@pattern))){
      blend <- 5 # blendPlus
    }
  } else{
    values <- c(0)
  }
  struct@pattern[struct@pattern==0] <- NA

  temp <- morphCpp(vals = vals,
                   kernel = struct@pattern,
                   value = values,
                   valRows = dims[2],
                   valCols = dims[1],
                   blend = blend,
                   merge = 2,
                   rotateKernel = FALSE,
                   strictKernel = FALSE)

  # build output ----
  .update_mosaik(obj, values = temp, step = step)
}
