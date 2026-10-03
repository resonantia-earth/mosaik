#' Offset the values in a mosaik
#'
#' Apply arithmetic operators to offset all cell values.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param fun [`character(1)`][character]\cr function used for offsetting; supported
#'   are all \link[=^]{arithmetic operators}.
#' @param value [`numeric(1)`][numeric]\cr the value by which to offset all
#'   values (can also be negative).
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}.
#' @examples
#' # the canopy height plus 5 and times 2, and the land cover modulo 3
#' m <- landscape |>
#'   mdf_offset(fun = "+", value = 5, layer = "canopy", add = "plus_5") |>
#'   mdf_offset(fun = "*", value = 2, layer = "canopy", add = "times_2") |>
#'   mdf_offset(fun = "%%", value = 3, layer = "cover", add = "modulo_3")
#' msk_vis(m, .layer("canopy"), .layer("plus_5"), .layer("times_2"))
#' msk_vis(m, .layer("cover"), .layer("modulo_3"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertChoice assertNumber assertCharacter
#' @export

mdf_offset <- function(obj = NULL,
                       fun = "+",
                       value = 1,
                       layer = NULL,
                       add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumber(x = value, finite = TRUE)
  assertChoice(x = fun, choices = c("+", "-", "*", "/", "%%", "%/%", "^", "**"))
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  temp <- do.call(what = fun, args = list(vals, value))

  # build output ----
  .update_mosaik(obj, values = temp, step = step)
}
