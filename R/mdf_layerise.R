#' Segregate values of a mosaik into layers
#'
#' Each distinct value in the source layer becomes its own layer in the output
#' mosaik. Cells belonging to a value are set to 1 (or their original value
#' when \code{flatten = FALSE}); all other cells in that layer are set to
#' \code{background}. Layers are stored with RLE compression, so the space
#' overhead is modest even for many classes.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param by [character(1)][character]\cr optional name of a layer in \code{obj}
#'   whose distinct values define the split. If left empty, the values of
#'   \code{layer} are used. To split by a layer of another mosaik, bring it in
#'   with \code{\link{msk_add}} first.
#' @param flatten [logical(1)][logical]\cr whether all values should be set to 1
#'   or the original value should be retained.
#' @param background [integerish(1)][integer]\cr the value any cell with value
#'   NA should have.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @return A mosaik with one layer per unique value, named after the values.
#' @examples
#' mdf_layerise(landscape, layer = "cover")
#' mdf_layerise(landscape, flatten = TRUE, layer = "cover")
#' @family operators to modify the overall object
#' @importFrom checkmate assertClass assertLogical assertIntegerish
#'   assertCharacter
#' @export

mdf_layerise <- function(obj = NULL,
                         by = NULL,
                         flatten = FALSE,
                         background = NA,
                         layer = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertLogical(x = flatten, len = 1)
  assertIntegerish(x = background)
  assertCharacter(x = layer, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # if 'by' is given, use its values for subsetting
  if(!is.null(by)){
    subVals <- .pull_layer(obj, by, "by")
    uVals <- sort(unique(subVals[!is.na(subVals)]))
  } else {
    subVals <- vals
    uVals <- sort(unique(vals[!is.na(vals)]))
  }

  # body ----
  new_layers <- list()
  for(i in seq_along(uVals)){
    temp <- rep(background, length(vals))
    mask <- subVals == uVals[i] & !is.na(subVals)
    if(flatten){
      temp[mask] <- 1L
    } else {
      temp[mask] <- vals[mask]
    }
    new_layers[[as.character(uVals[i])]] <- temp
  }

  # build output ----
  prov <- msk_prov("mdf_layerise", list(flatten = flatten,
                     background = background, layer = layer))

  new_mosaik(
    extent     = obj@extent,
    dims       = obj@dims,
    layers     = new_layers,
    categories = list(),
    patches    = list(),
    global     = list(),
    crs        = obj@crs,
    provenance = c(obj@provenance, list(prov))
  )
}
