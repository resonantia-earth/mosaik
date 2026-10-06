#' Blend layers of a mosaik
#'
#' Combine two or more layers of a mosaik into one, cell by cell.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param layers [`character(.)`][character]\cr the layers to blend, at least two.
#'   Defaults to all layers of \code{obj}, in the order they are stored.
#' @param fun [`character(1)`][character] or [`function`][function]\cr how the
#'   values of a cell are combined. Either an arithmetic operator written as
#'   text (\code{"+"}, \code{"-"}, \code{"*"}, \code{"/"}, \code{"\%\%"},
#'   \code{"\%/\%"}, \code{"^"}), applied from the first layer to the last; one
#'   of the summaries shared with \code{\link{mdf_summarise}}: \code{"max"},
#'   \code{"min"}, \code{"sum"}, \code{"mean"}, \code{"median"}, \code{"any"}
#'   (any non-zero value present), \code{"all"}, \code{"n"} (number of layers),
#'   \code{"n_distinct"} (number of distinct non-zero values) or
#'   \code{"unique"} (the one non-zero value present, \code{NA} if there are
#'   several); or a function that takes the values of one cell, one per layer,
#'   and returns a single value.
#' @param weights [`numeric(.)`][numeric]\cr optional weights, one per layer, by
#'   which values are multiplied before \code{fun} is applied.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
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
#'   An operator works on whole layers at once and is fast. A summary or a
#'   function is evaluated once per cell, so it is slower.
#' @examples
#' # the distance to the forest, to the river, and to whichever is nearer
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_filter(cover == 1, add = "river") |>
#'   mdf_distance(layer = "forest", add = "to_forest") |>
#'   mdf_distance(layer = "river", add = "to_river") |>
#'   mdf_blend(layers = c("to_forest", "to_river"), fun = "min",
#'             add = "to_either")
#' msk_vis(m, .layer("to_forest"), .layer("to_river"), .layer("to_either"))
#'
#' # the forest as 1, the river as 2, both counted: a weighted sum of 0/1 layers
#' m <- mdf_blend(m, layers = c("forest", "river"), fun = "+",
#'                weights = c(1, 2), add = "forest_river")
#' msk_vis(m, .layer("forest"), .layer("river"), .layer("forest_river"))
#'
#' # a layer held in another mosaik is brought in first: here the canopy height
#' # is damped by a gradient that falls from west to east
#' other <- mosaik(extent = c(0, 60, 0, 56), res = 1) |>
#'   drw_gradient(type = "planar", invert = TRUE, name = "damping")
#' m <- msk_add(landscape, other, damping) |>
#'   mdf_blend(layers = c("canopy", "damping"), fun = "*", add = "damped")
#' msk_vis(m, .layer("canopy"), .layer("damping"), .layer("damped"))
#' @family operators to modify the overall object
#' @importFrom checkmate assertClass assertCharacter assertNumeric
#' @export

mdf_blend <- function(obj = NULL,
                      layers = NULL,
                      fun = "+",
                      weights = NULL,
                      add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

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
  operators <- c("+", "-", "*", "/", "%%", "%/%", "^")
  if(is.character(fun) && length(fun) == 1 && fun %in% operators){
    temp <- allVals[[1]]
    for(i in 2:length(allVals)){
      temp <- do.call(what = fun, args = list(temp, allVals[[i]]))
    }
  } else {
    f <- .summary_fun(fun, "mdf_blend")
    temp <- as.numeric(reduceCpp(lVals = allVals, f = f))
  }

  # build output ----
  .update_mosaik(obj, values = temp, keep = FALSE, step = step)
}
