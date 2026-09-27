#' Blend layers of a mosaik
#'
#' Combine two or more layers of a mosaik into one, cell by cell.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param layers [character(.)][character]\cr the layers to blend, at least two.
#'   Defaults to all layers of \code{obj}, in the order they are stored.
#' @param fun [character(1) or function(1)][character]\cr how values are
#'   combined. For character, one of the arithmetic operators
#'   (\code{"+", "-", "*", "/", "\%\%", "\%/\%", "^", "**"}), applied
#'   successively from left to right. For function, a function taking a numeric
#'   vector of one value per layer and returning a scalar (evaluated per cell
#'   via \code{reduceCpp}).
#' @param weights [numeric(.)][numeric]\cr optional weights, one per layer, by
#'   which values are multiplied before \code{fun} is applied.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   the first blended layer; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}, carrying the blended
#'   layer.
#' @details Blending happens within one object, so layers computed elsewhere are
#'   brought in with \code{\link{msk_add}} first. That keeps the grid check in
#'   one place and makes the inputs of a blend visible as named layers:
#'
#'   \preformatted{
#'   obj |>
#'     msk_add(other, "suitability") |>
#'     mdf_blend(layers = c("cover", "suitability"), fun = "*")
#'   }
#'
#'   Character \code{fun} is the fast path and covers most cases. A function is
#'   for reductions no operator expresses, such as \code{max} or \code{var}
#'   across layers.
#' @examples
#' # blend two named layers
#' mdf_blend(landscape, layers = c("cover", "intensity"), fun = "*")
#'
#' # blend every layer, weighted, into a new layer
#' mdf_blend(landscape, fun = "+", weights = c(0.3, 0.7), add = "score")
#'
#' # a reduction that no arithmetic operator expresses
#' mdf_blend(landscape, fun = max)
#'
#' # a layer held in another mosaik is brought in first, then blended
#' other <- mosaik(extent = c(0, 60, 0, 56), res = 1,
#'                 vals = list(elevation = runif(60 * 56, 0, 800)))
#' msk_add(landscape, other, elevation) |>
#'   mdf_blend(layers = c("intensity", "elevation"), fun = "+", add = "sum")
#' @family operators to modify the overall object
#' @importFrom checkmate assertClass assertCharacter assertNumeric
#' @export

mdf_blend <- function(obj = NULL,
                      layers = NULL,
                      fun = "+",
                      weights = NULL,
                      add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layers, min.len = 2, null.ok = TRUE, any.missing = FALSE)
  assertNumeric(x = weights, null.ok = TRUE, any.missing = FALSE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # resolve layers ----
  if(is.null(layers)) layers <- names(obj@layers)
  if(length(layers) < 2){
    stop("'obj' has only one layer; nothing to blend.", call. = FALSE)
  }
  missing <- setdiff(layers, names(obj@layers))
  if(length(missing) > 0){
    stop("layer(s) not found in 'obj': ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  # pull data ----
  allVals <- lapply(layers, function(nm) msk_pull(obj, nm))

  if(!is.null(weights)){
    if(length(weights) != length(allVals)){
      stop("'weights' must have one value per blended layer.", call. = FALSE)
    }
    allVals <- lapply(seq_along(allVals), function(i) allVals[[i]] * weights[i])
  }

  # body ----
  if(is.character(fun)){
    temp <- allVals[[1]]
    for(i in 2:length(allVals)){
      temp <- do.call(what = fun, args = list(temp, allVals[[i]]))
    }
  } else if(is.function(fun)){
    temp <- as.numeric(reduceCpp(lVals = allVals, f = fun))
  } else {
    stop("'fun' must be a character (arithmetic operator) or a function.",
         call. = FALSE)
  }

  # build output ----
  out_layer <- .resolve_add(obj, layers[1], add)
  prov <- msk_prov("mdf_blend",
                     list(fun = if(is.character(fun)) fun else "custom",
                          inputs = paste(layers, collapse = ", "),
                          layer = out_layer))
  msk_set(obj, out_layer, temp, prov, keep = FALSE)
}
