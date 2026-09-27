#' Select cells by a predicate over layers
#'
#' Evaluate a logical expression over the layers of a mosaik and write the
#' result back into a layer, either keeping the cells that satisfy it or
#' recording the predicate itself as a mask.
#'
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param expr an expression evaluated in the context of layer values (e.g.
#'   \code{cover > 3}). Any layer of \code{obj} may be named.
#' @param value [logical(1)][logical] or [numeric(1)][numeric]\cr what to write where the predicate
#'   holds. \code{FALSE} (default) keeps \code{layer}'s own values and sets
#'   everything else to \code{background} -- a filter. \code{TRUE} writes the
#'   predicate itself as a 1/0 mask. A number writes that number instead of 1.
#' @param background [numeric(1)][numeric]\cr the value for cells where the
#'   predicate does not hold. Defaults to \code{NA}; pass \code{0} for a
#'   ready-made mask that needs no following \code{\link{mdf_replace}}.
#' @param layer [character(1)][character]\cr the layer to operate on. Defaults
#'   to the first layer. When \code{value} is set the layer is only the
#'   destination and its own values are not read.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}.
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
#'   Set operations want \code{value = TRUE, background = 0}, which yields a
#'   clean 1/0 mask. The default (\code{value = FALSE}) instead keeps the
#'   layer's values, which is what "filter this layer" means but is rarely what
#'   a set operation wants -- with the default, a predicate over other layers
#'   returns \code{layer}'s values at those cells, not the predicate.
#' @examples
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#'
#' # keep the values of cells that satisfy the predicate
#' mdf_filter(landscape, cover == 47)
#'
#' # record a set difference as a mask in a new layer
#' forest |>
#'   mdf_erode(add = "core") |>
#'   mdf_filter(cover == 1 & core == 0, value = TRUE, background = 0,
#'              add = "rim")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertCharacter assertNumber
#' @export

mdf_filter <- function(obj = NULL,
                       expr,
                       value = FALSE,
                       background = NA,
                       layer = NULL,
                       add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  if(is.null(layer)) layer <- names(obj@layers)[1]
  assertCharacter(x = layer, len = 1, null.ok = FALSE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)
  assertNumber(x = background, na.ok = TRUE)

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

  if(isFALSE(value)){
    vals <- msk_pull(obj, layer)
    vals[!hit] <- background
  } else {
    fill <- if(isTRUE(value)) 1 else value
    vals <- rep(background, prod(obj@dims))
    vals[hit] <- fill
  }

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_filter", list(expr = deparse(predicate),
                                        layer = out_layer))
  msk_set(obj, out_layer, vals, prov)
}
