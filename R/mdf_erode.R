#' Morphologically erode a mosaik
#'
#' The morphological operation 'erode' decreases the value of cells that match a
#' structuring element.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param struct [`struct(1)`][struct]\cr the structuring element; see
#'   \code{\link{msk_struct}} for details.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, where an erosion has
#'   been performed.
#' @examples
#' # the forest shrunk by one cell, with the default disc and a 3 x 3 square
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_erode(layer = "forest", add = "eroded") |>
#'   mdf_erode(struct = msk_struct("square", width = 3), layer = "forest",
#'             add = "eroded_square")
#' msk_vis(m, .layer("forest"), .layer("eroded"), .layer("eroded_square"))
#' @family operators to morphologically modify a raster
#' @importFrom checkmate assertClass assertCharacter
#' @export

mdf_erode <- function(obj = NULL,
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
      blend <- 6 # blendMinus
    }
  } else{
    values <- c(1)
  }
  struct@pattern[struct@pattern==0] <- NA

  temp <- morphCpp(vals = vals,
                   kernel = struct@pattern,
                   value = values,
                   valRows = dims[2],
                   valCols = dims[1],
                   blend = blend,
                   merge = 1,
                   rotateKernel = FALSE,
                   strictKernel = FALSE)

  # build output ----
  .update_mosaik(obj, values = temp, step = step)
}
