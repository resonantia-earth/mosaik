#' Morphologically modify a mosaik
#'
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param struct [`struct(1)`][struct]\cr the structuring element; see
#'   \code{\link{msk_struct}} for details.
#' @param blend [`character(1)`][character]\cr \code{identity}, \code{equal},
#'   \code{lower}, \code{greater}, \code{plus}, \code{minus}, \code{product};
#'   see Details.
#' @param merge [`character(1)`][character]\cr \code{min}, \code{max}, \code{all},
#'   \code{any}, \code{!all}, \code{!any}, \code{sum}, \code{mean},
#'   \code{median}, \code{sd}, \code{cv}, \code{sumNa}; see Details.
#' @param strict [`logical(1)`][logical]\cr how cells at the edge of the grid are
#'   treated, where the kernel reaches beyond it. \code{FALSE} (default): the
#'   kernel is clipped to the grid, so the operation covers every cell.
#'   \code{TRUE}: only cells whose kernel fits entirely inside the grid are
#'   computed; the others get \code{background}.
#' @param rotate [`logical(1)`][logical]\cr whether to try all four quarter turns
#'   of the kernel and stop at the first whose result is 1. This is meant for
#'   matching a pattern in any orientation (see \code{\link{mdf_match}}), and
#'   needs a square kernel. Default \code{FALSE}.
#' @param background [`integerish(1)`][integer]\cr the value any cell with value
#'   NA should have.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @details The \code{morphCpp} function (internal) iteratively goes through
#'   each pixel and compares a structuring element with the mosaik at that
#'   location. The structuring element is a small grid of values centred on the
#'   cell: \code{NA} cells are outside the neighbourhood. A \code{0} means
#'   "this neighbour must be 0" for \code{blend = "equal"} (matching a
#'   pattern) and is outside the neighbourhood for every other blend, so a
#'   disc acts as a disc. The result depends on \code{blend} and \code{merge}:
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
#' # smoothing (the mean of the neighbourhood) and the local maximum
#' m <- landscape |>
#'   mdf_morph(struct = s, blend = "product", merge = "mean",
#'             layer = "canopy", add = "mean") |>
#'   mdf_morph(struct = s, blend = "identity", merge = "max",
#'             layer = "canopy", add = "max")
#' msk_vis(m, .layer("canopy"), .layer("mean"), .layer("max"))
#'
#' # the interior of a binary layer: cells whose whole neighbourhood is 1
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_morph(struct = s, blend = "equal", merge = "all", layer = "forest",
#'             add = "interior")
#' msk_vis(m, .layer("forest"), .layer("interior"))
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
                      rotate = FALSE,
                      strict = FALSE,
                      background = NA,
                      layer = NULL,
                      add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

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
  if (rotate && nrow(struct@pattern) != ncol(struct@pattern)) {
    stop("'rotate = TRUE' needs a square 'struct'; this one is ",
         nrow(struct@pattern), " by ", ncol(struct@pattern), " cells.",
         call. = FALSE)
  }
  assertIntegerish(x = background)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  uVals <- unique(vals)

  # body ----
  # a 0 in the pattern means "must be 0" only when matching (equal); for any
  # other blend it marks a cell outside the neighbourhood, as in mdf_dilate,
  # so a disc acts as a disc and not as the square around it
  kernel <- struct@pattern
  if (blend != "equal") kernel[kernel == 0] <- NA

  temp <- morphCpp(vals = vals,
                   valRows = dims[2],
                   valCols = dims[1],
                   kernel = kernel,
                   value = uVals,
                   blend = blendID,
                   merge = mergeID,
                   rotateKernel = rotate,
                   strictKernel = strict)
  temp[is.na(temp)] <- background

  # build output ----
  .update_mosaik(obj, values = temp, step = step)
}
