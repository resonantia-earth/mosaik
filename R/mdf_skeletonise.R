#' Determine the skeleton of foreground patches
#'
#' The morphological skeleton preserves the extent and connectivity (topology)
#' of the patch.
#' @param obj [mosaik]\cr a binarised mosaik to modify.
#' @param background [integerish(1)][integer]\cr the value any cell with value
#'   NA should have.
#' @param anchor [character(1)][character]\cr the name of a layer whose non-zero
#'   cells anchor the skeleton: they are never removed during thinning, so the
#'   skeleton stays connected to them and its free ends rest on them. If
#'   \code{NULL} (default) the ordinary unanchored skeleton is computed.
#' @param method [character(1)][character]\cr the thinning algorithm,
#'   \code{"zhangSuen"} (default) or \code{"homotopic"}.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}, in which foreground
#'   patches have been transformed into their morphological skeletons.
#' @details \code{"zhangSuen"} thins in two directional sub-passes, so which
#'   cell of an even-width strand survives depends on the scan direction.
#'   \code{"homotopic"} instead deletes simple points -- cells whose removal
#'   keeps foreground and background connected as before -- independently of
#'   that order, and is what MSPA builds on (see
#'   \code{vignette("mspa", package = "mosaik")}).
#' @examples
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' mdf_skeletonise(forest)
#' mdf_skeletonise(forest, method = "homotopic")
#' @family operators to determine objects
#' @references Ranwez, V. & Soille, P. (2002). Order independent homotopic
#'   thinning for binary and grey tone anchored skeletons. \emph{Pattern
#'   Recognition Letters} 23(6), 687-702.
#' @importFrom checkmate assertClass assertIntegerish assertCharacter
#'   assertChoice
#' @export

mdf_skeletonise <- function(obj = NULL,
                            background = NA,
                            anchor = NULL,
                            method = "zhangSuen",
                            layer = NULL,
                            add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertIntegerish(x = background)
  assertCharacter(x = anchor, len = 1, null.ok = TRUE)
  assertChoice(x = method, choices = c("zhangSuen", "homotopic"))
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims

  # body ----
  if(!isBinaryCpp(vals = vals)){
    stop("'obj' is not binary, please run 'mdf_binarise()' first.")
  }
  anchorVec <- if(is.null(anchor)) numeric(0) else as.numeric(msk_pull(obj, anchor))
  temp <- skeletoniseCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                         anchor = anchorVec,
                         homotopic = method == "homotopic")
  temp[temp == 0] <- background

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_skeletonise", list(background = background, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
