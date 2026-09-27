#' Morphologically modify a mosaik
#'
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param struct [struct(1)][struct]\cr the structuring element; see
#'   \code{\link{msk_struct}} for details.
#' @param blend [character(1)][character]\cr \code{identity}, \code{equal},
#'   \code{lower}, \code{greater}, \code{plus}, \code{minus}, \code{product};
#'   see Details.
#' @param merge [character(1)][character]\cr \code{min}, \code{max}, \code{all},
#'   \code{any}, \code{sum}, \code{mean}, \code{median}, \code{sd}, \code{cv},
#'   \code{one}, \code{zero}, \code{na}; see Details.
#' @param strict [logical(1)][logical]\cr whether ...
#' @param rotate [logical(1)][logical]\cr whether the kernel should be rotated.
#' @param background [integerish(1)][integer]\cr the value any cell with value
#'   NA should have.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @details The \code{morphCpp} function (internal) iteratively goes through
#'   each pixel and compares a structuring element with the mosaik at that
#'   location. The result depends on \code{blend} and \code{merge}:
#'   \itemize{
#'     \item First, values covered by the structuring element are blended
#'       pairwise with the element (\code{blend}).
#'     \item Then these values are merged into a single value (\code{merge})
#'       and assigned in the current location.
#'   }
#'
#'   Blend functions: identity, equal, lower, greater, plus, minus, product.
#'   Merge functions: min, max, all, any, !all, !any, sum, mean, median, sd,
#'   cv, sumNa.
#' @examples
#' s <- msk_struct("disc", width = 3)
#'
#' # smoothing: weighted mean of neighbourhood
#' mdf_morph(landscape, struct = s, blend = "product", merge = "mean",
#'           layer = "intensity")
#'
#' # edge detection on binary layer
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' mdf_morph(forest, struct = s, blend = "equal", merge = "all",
#'           strict = TRUE, background = 0)
#'
#' # local maximum
#' mdf_morph(landscape, struct = s, blend = "identity", merge = "max",
#'           layer = "intensity")
#' @references Credit for the original idea and C++ code is due to Jon Clayden
#'   (\href{https://github.com/jonclayden/mmand}{R::mmand}).
#' @family operators to morphologically modify a raster
#' @importFrom checkmate assertClass assertSubset assertLogical assertIntegerish
#'   assertCharacter
#' @export

mdf_morph <- function(obj = NULL,
                      struct = NULL,
                      blend = NULL,
                      merge = NULL,
                      rotate = TRUE,
                      strict = TRUE,
                      background = NA,
                      layer = NULL,
                      add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertClass(x = struct, classes = "struct", null.ok = TRUE)
  if(is.null(struct)) struct <- msk_struct(type = "diamond", width = 3, height = 3)
  assertSubset(blend, choices = c("identity", "equal", "lower", "greater", "plus", "minus", "product"), empty.ok = FALSE)
  blendID <- which(c("identity", "equal", "lower", "greater", "plus", "minus", "product") %in% blend)
  assertSubset(merge, choices = c("min", "max", "all", "any", "!all", "!any", "sum", "mean", "median", "sd", "cv", "sumNa"), empty.ok = FALSE)
  mergeID <- which(c("min", "max", "all", "any", "!all", "!any", "sum", "mean", "median", "sd", "cv", "sumNa") %in% merge)
  assertLogical(x = rotate, any.missing = FALSE)
  assertLogical(x = strict, any.missing = FALSE)
  assertIntegerish(x = background)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  uVals <- unique(vals)

  # body ----
  temp <- morphCpp(vals = vals,
                   valRows = dims[2],
                   valCols = dims[1],
                   kernel = struct@pattern,
                   value = uVals,
                   blend = blendID,
                   merge = mergeID,
                   rotateKernel = rotate,
                   strictKernel = strict)
  temp[is.na(temp)] <- background

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_morph", list(struct = struct@pattern, blend = blend,
                     merge = merge, rotate = rotate, strict = strict,
                     background = background, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
