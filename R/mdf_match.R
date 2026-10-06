#' Select cells based on a structuring element
#'
#' Transform a mosaik by setting all cells that do not match a structuring
#' element to \code{background}. Also known as the 'hit-or-miss'-transform.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param struct [`struct(1)`][struct]\cr the pattern to match; see
#'   \code{\link{msk_struct}} for details.
#' @param rotate [`logical(1)`][logical]\cr whether the kernel should be rotated
#'   for all possible rotations or used as is.
#' @param background [`integerish(1)`][integer]\cr the value any non-matching cell
#'   should have.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}.
#' @examples
#' # a forest cell with no forest around it, and the end of a forest line one
#' # cell wide: the line leaves the end downwards, the cells marked NA may be
#' # anything
#' isolated <- msk_struct(custom = matrix(c(0, 0, 0,
#'                                          0, 1, 0,
#'                                          0, 0, 0), nrow = 3, byrow = TRUE))
#' end <- msk_struct(custom = matrix(c( 0, 0,  0,
#'                                      0, 1,  0,
#'                                     NA, 1, NA), nrow = 3, byrow = TRUE))
#'
#' # isolated cells, line ends in all four directions, and only those whose line
#' # leaves downwards; drawn over the forest
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_match(struct = isolated, layer = "forest", add = "isolated") |>
#'   mdf_match(struct = end, layer = "forest", add = "ends") |>
#'   mdf_match(struct = end, rotate = FALSE, layer = "forest", add = "ends_down")
#' over <- function(layer) {
#'   list(.layer("forest", panel = layer, colours = c("grey92", "grey70"),
#'               legend = FALSE),
#'        .layer(layer, panel = layer, colours = c("black", "black"),
#'               legend = FALSE))
#' }
#' do.call(msk_vis, c(list(m), over("isolated"), over("ends"), over("ends_down")))
#' @family operators to select a subset of cells
#' @export

mdf_match <- function(obj = NULL,
                      struct = NULL,
                      rotate = TRUE,
                      background = NA,
                      layer = NULL,
                      add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  out <- mdf_morph(obj = obj,
                   struct = struct,
                   blend = "equal",
                   merge = "sumNa",
                   rotate = rotate,
                   strict = TRUE,
                   background = background,
                   layer = layer,
                   add = add)

  .update_mosaik(out, step = step)
}
