#' Morphologically erode a mosaik
#'
#' The morphological operation 'erode' decreases the value of cells that match a
#' structuring element.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param struct [struct(1)][struct]\cr the structuring element; see
#'   \code{\link{msk_struct}} for details.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, where an erosion has
#'   been performed.
#' @examples
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' mdf_erode(forest)
#' mdf_erode(forest, struct = msk_struct("square", width = 3))
#' @family operators to morphologically modify a raster
#' @importFrom checkmate assertClass assertCharacter
#' @export

mdf_erode <- function(obj = NULL,
                      struct = NULL,
                      layer = NULL,
                      add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

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
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_erode", list(struct = struct@pattern, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
