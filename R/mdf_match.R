#' Select cells based on a structuring element
#'
#' Transform a mosaik by setting all cells that do not match a structuring
#' element to \code{background}. Also known as the 'hit-or-miss'-transform.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param struct [struct(1)][struct]\cr the pattern to match; see
#'   \code{\link{msk_struct}} for details.
#' @param rotate [logical(1)][logical]\cr whether the kernel should be rotated
#'   for all possible rotations or used as is.
#' @param background [integerish(1)][integer]\cr the value any non-matching cell
#'   should have.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}.
#' @examples
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' s <- msk_struct("cross", width = 3)
#' mdf_match(forest, struct = s)
#' mdf_match(forest, struct = s, rotate = FALSE)
#' @family operators to select a subset of cells
#' @export

mdf_match <- function(obj = NULL,
                      struct = NULL,
                      rotate = TRUE,
                      background = NA,
                      layer = NULL,
                      add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  out <- mdf_morph(obj = obj,
                   struct = struct,
                   blend = "equal",
                   merge = "sumNa",
                   rotate = rotate,
                   background = background,
                   layer = layer,
                   add = add)

  return(out)
}
