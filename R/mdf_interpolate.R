#' Spatially smooth or filter cell values
#'
#' Apply a weighted spatial filter (convolution) to a layer. Each cell is
#' replaced by the weighted mean of its neighbourhood defined by \code{struct}.
#' This is a convenience wrapper around \code{\link{mdf_morph}} with
#' \code{blend = "product"} and \code{merge = "mean"}.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param struct [struct(1)][struct]\cr the structuring element whose values
#'   serve as weights. Defaults to a 3x3 disc (uniform weights). See
#'   \code{\link{msk_struct}}.
#' @param layer [character(1)][character]\cr the layer to filter. Defaults to
#'   the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions with smoothed cell values.
#' @examples
#' mdf_interpolate(landscape, layer = "intensity")
#' mdf_interpolate(landscape, struct = msk_struct("square", width = 5),
#'                 layer = "intensity")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertCharacter
#' @export

mdf_interpolate <- function(obj = NULL,
                            struct = NULL,
                            layer = NULL,
                            add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  if(is.null(struct)) struct <- msk_struct(type = "diamond", width = 3, height = 3)

  # delegate to mdf_morph ----
  out <- mdf_morph(obj = obj,
                   struct = struct,
                   blend = "product",
                   merge = "mean",
                   rotate = FALSE,
                   strict = FALSE,
                   background = NA,
                   layer = layer,
                   add = add)

  # fix provenance to reflect mdf_interpolate
  out@provenance[[length(out@provenance)]]$fn <- "mdf_interpolate"

  out
}
