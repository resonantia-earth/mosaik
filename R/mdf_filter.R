#' Select cells by a predicate over layers
#'
#' Evaluate a logical expression over the layers of a mosaik and write the
#' result into a layer, either as a mask of the cells that satisfy it or as the
#' values of \code{layer} at those cells.
#'
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param expr an expression evaluated in the context of layer values (e.g.
#'   \code{cover > 3}). Any layer of \code{obj} may be named.
#' @param value [`logical(1)`][logical]\cr what to write. \code{FALSE} (default)
#'   writes the predicate as a mask: 1 where it holds, 0 elsewhere. \code{TRUE}
#'   keeps the values of \code{layer} where the predicate holds and writes
#'   \code{NA} elsewhere.
#' @param layer [`character(1)`][character]\cr the layer whose values are kept
#'   when \code{value = TRUE}. Defaults to the first layer. With
#'   \code{value = FALSE} it is only the destination when \code{add} is
#'   \code{NULL}.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik with the mask or the kept values.
#' @details
#'   Because the predicate may name any layer, this one operator is the whole
#'   boolean algebra over layers, and mosaik therefore has no separate union,
#'   intersection or difference operators:
#'
#'   \describe{
#'     \item{difference}{\code{a == 1 & b == 0}}
#'     \item{union}{\code{a == 1 | b == 1}}
#'     \item{intersection}{\code{a == 1 & b == 1}}
#'     \item{complement}{\code{a == 0}}
#'   }
#'
#'   Cells where the predicate is \code{NA} count as not satisfying it.
#' @examples
#' # the forest as a mask, and the canopy height of the forest cells only
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_filter(cover == 47, value = TRUE, layer = "canopy",
#'              add = "forest_canopy")
#' msk_vis(m, .layer("cover"), .layer("forest"), .layer("forest_canopy"))
#'
#' # a set difference: the forest cells that are not core, i.e. its rim
#' m <- m |>
#'   mdf_erode(layer = "forest", add = "core") |>
#'   mdf_filter(forest == 1 & core == 0, add = "rim")
#' msk_vis(m, .layer("forest"), .layer("core"), .layer("rim"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertCharacter assertFlag
#' @export

mdf_filter <- function(obj = NULL,
                       expr,
                       value = FALSE,
                       layer = NULL,
                       add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertFlag(x = value)
  if(is.null(layer)) layer <- names(obj@layers)[1]
  assertCharacter(x = layer, len = 1, null.ok = FALSE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  env <- lapply(names(obj@layers), function(nm) msk_pull(obj, nm))
  names(env) <- names(obj@layers)

  # body ----
  # Normally the predicate arrives unevaluated and substitute() captures it.
  # When a recipe is replayed, do.call passes the stored language object by
  # value, so substitute() returns a promise-free symbol bound to it; unwrap
  # that once so both paths end up with the same expression.
  predicate <- substitute(expr)
  if (is.name(predicate) &&
      exists(as.character(predicate), envir = parent.frame(), inherits = FALSE)) {
    bound <- get(as.character(predicate), envir = parent.frame(), inherits = FALSE)
    if (is.call(bound) || is.name(bound)) predicate <- bound
  }
  mask <- eval(predicate, envir = env, enclos = parent.frame())
  if(!is.logical(mask) || length(mask) != prod(obj@dims)){
    stop("predicate must evaluate to a logical vector of length ", prod(obj@dims))
  }
  hit <- mask & !is.na(mask)

  if(value){
    vals <- msk_pull(obj, layer)
    vals[!hit] <- NA
  } else {
    vals <- as.integer(hit)
  }

  # build output ----
  # a mask is a new kind of value; kept values stay what they were
  .update_mosaik(obj, values = vals, keep = value, step = step)
}
