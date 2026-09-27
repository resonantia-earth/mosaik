#' Offset the values in a mosaik
#'
#' Apply arithmetic operators to offset all cell values.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param fun [character(1)][character]\cr function used for offsetting; supported
#'   are all \link[=^]{arithmetic operators}.
#' @param value [numeric(1)][numeric]\cr the value by which to offset all
#'   values (can also be negative).
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}.
#' @examples
#' mdf_offset(landscape, fun = "+", value = 10, layer = "intensity")
#' mdf_offset(landscape, fun = "*", value = 2, layer = "intensity")
#' mdf_offset(landscape, fun = "%%", value = 3, layer = "cover")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertChoice assertIntegerish
#'   assertCharacter
#' @export

mdf_offset <- function(obj = NULL,
                       fun = "+",
                       value = 1,
                       layer = NULL,
                       add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertIntegerish(x = value, any.missing = FALSE, len = 1)
  assertChoice(x = fun, choices = c("+", "-", "*", "/", "%%", "%/%", "^", "**"))
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  temp <- do.call(what = fun, args = list(vals, value))

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_offset", list(fun = fun, value = value, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
