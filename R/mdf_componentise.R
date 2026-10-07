#' Determine patches in a mosaik
#'
#' A patch is a set of connected cells that share a value. Each patch gets its
#' own number; cells with the value 0 or \code{NA} belong to no patch and get
#' \code{background}. On a binary layer, the patches are those of the
#' foreground (1); on a categorical layer, such as land cover, those of every
#' class.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param connectivity [`integer(1)`][integer]\cr neighbourhood connectivity;
#'   \code{4} for rook/diamond (default) or \code{8} for queen/box.
#' @param background [`integerish(1)`][integer]\cr the value of every cell that
#'   is not in a patch, i.e. cells with the value 0 or \code{NA}. Keep the
#'   default \code{NA} to measure the patches: any other value is a class of
#'   the layer and is measured as one more group.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik in which each cell carries the number of its patch.
#' @details The layer of patch numbers is a grouping: each of its classes is a
#'   patch, so the \code{msr_*} functions and \code{\link{msr}} measure the
#'   patches when given this layer. How cells connect into patches is set here
#'   and never decided by a measure. Which class a patch came from is read
#'   from the source layer through the patch's cells, for example
#'   \code{msr(m, "mean(cover.self == 47)", label = "forest", layer = "patch")}.
#' @examples
#' # the forest patches with rook (4) and queen (8) neighbours; drawn in shuffled
#' # colours, since neighbouring patch numbers would otherwise get similar ones
#' set.seed(1)
#' shuffled <- function(m, layer) {
#'   sample(hcl.colors(max(msk_pull(m, layer), na.rm = TRUE), "Dark 3"))
#' }
#' f <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_componentise(layer = "forest", add = "rook") |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "queen")
#' msk_vis(f, .layer("forest"),
#'         .layer("rook", colours = shuffled(f, "rook"), legend = FALSE),
#'         .layer("queen", colours = shuffled(f, "queen"), legend = FALSE))
#'
#' # the patches of every land-cover class, and the class of each patch
#' m <- mdf_componentise(landscape, layer = "cover", add = "patch") |>
#'   msr(equation = "cover.self[1]", label = "source", layer = "patch")
#' table(msk_table(m, layer = "patch")$source)
#' msk_vis(m, .layer("cover"),
#'         .layer("patch", colours = shuffled(m, "patch"), legend = FALSE))
#' @family operators to determine objects
#' @importFrom checkmate assertClass assertIntegerish assertChoice
#'   assertCharacter
#' @export

mdf_componentise <- function(obj = NULL,
                             connectivity = 4L,
                             background = NA,
                             layer = NULL,
                             add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

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
  # connected cells of equal value form a patch; 0 and NA form none
  fg <- vals
  fg[!is.na(fg) & fg == 0] <- NA
  temp <- componentsCpp(vals = fg, nrow = dims[2], ncol = dims[1],
                        connectivity = connectivity)
  temp[is.na(temp)] <- background

  # build output ----
  # componentising changes the layer's kind (-> patch IDs); drop any prior role
  .update_mosaik(obj, values = temp, keep = FALSE, step = step)
}
