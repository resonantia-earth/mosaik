#' Select cells based on a mask
#'
#' Transform a mosaik by setting all cells that are not covered by a mask to
#' \code{background}.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param by [character(1)][character]\cr name of a binary layer in \code{obj}
#'   where cells to retain have value 1 and all others value 0. To mask by a
#'   layer held in another mosaik, bring it in with \code{\link{msk_add}} first.
#' @param background [integerish(1)][integer]\cr the value any masked cell
#'   should have (default \code{NA}).
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}, where masked cells
#'   have been set to \code{background}.
#' @examples
#' # mask one layer by a binary layer of the same mosaik
#' m <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest")
#' mdf_mask(m, by = "forest", layer = "intensity")
#' mdf_mask(m, by = "forest", layer = "intensity", background = 0)
#' @family operators to select a subset of cells
#' @importFrom checkmate assertClass assertIntegerish assertCharacter
#' @export

mdf_mask <- function(obj = NULL,
                     by,
                     background = NA,
                     layer = NULL,
                     add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertIntegerish(x = background, null.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  maskVals <- .pull_layer(obj, by, "by")

  # body ----
  if(!isBinaryCpp(vals = maskVals)){
    stop("'by' is not binary, please run 'mdf_binarise()' first.")
  }

  temp <- vals
  temp[maskVals == 0] <- background

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_mask", list(background = background, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
