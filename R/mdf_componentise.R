#' Determine patches in a mosaik
#'
#' Patches are sets of cells that are connected.
#' @param obj [mosaik]\cr a binarised mosaik to modify.
#' @param connectivity [integer(1)][integer]\cr neighbourhood connectivity;
#'   \code{4} for rook/diamond (default) or \code{8} for queen/box.
#' @param background [integerish(1)][integer]\cr the value any cell with value
#'   NA should have.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, in which neighboring
#'   cells of the foreground have been assigned the same value, forming patches.
#' @examples
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' mdf_componentise(forest)
#' mdf_componentise(forest, connectivity = 8L)
#' @family operators to determine objects
#' @importFrom checkmate assertClass assertIntegerish assertChoice
#'   assertCharacter
#' @export

mdf_componentise <- function(obj = NULL,
                             connectivity = 4L,
                             background = NA,
                             layer = NULL,
                             add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = connectivity, choices = c(4L, 8L))
  assertIntegerish(x = background)
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
  temp <- componentsCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                        connectivity = connectivity)

  # set na
  temp[is.na(temp)] <- background

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_componentise", list(connectivity = connectivity,
                     background = background, layer = out_layer))
  # componentising changes the layer's kind (-> patch IDs); drop any prior role
  msk_set(obj, out_layer, temp, prov, keep_role = FALSE)
}
