#' Categorise a mosaik
#'
#' Transform the values of a mosaik to a (smaller) set of categories.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param breaks [numeric(.)][numeric]\cr the break points where categories
#'   should be delimited.
#' @param n [integerish(1)][integer]\cr number of categories.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, in which cells have
#'   the category value into which their values fall.
#' @details Using \code{n} will determine \code{breaks} based on the value
#'   range of \code{obj} so that values are assigned to n categories.
#' @examples
#' mdf_categorise(landscape, n = 3, layer = "intensity")
#' mdf_categorise(landscape, breaks = c(25, 50, 75), layer = "intensity")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertNumeric assertIntegerish
#'   assertCharacter
#' @export

mdf_categorise <- function(obj = NULL,
                           breaks = NULL,
                           n = NULL,
                           layer = NULL,
                           add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumeric(x = breaks, null.ok = TRUE)
  assertIntegerish(x = n, null.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  theRange <- range(vals, na.rm = TRUE)
  if(is.null(breaks)){
    assertIntegerish(x = n)
    breaks <- seq(theRange[1], theRange[2], length.out = n+1)
  }

  # manage the tails
  if(!any(breaks == theRange[1])){
    breaks <- c(theRange[1], breaks)
  }
  if(!any(breaks == theRange[2])){
    breaks <- c(breaks, theRange[2])
  }
  temp <- findInterval(vals, breaks, rightmost.closed = TRUE)

  # build output ----
  uVals <- sort(unique(temp[!is.na(temp)]))
  groupsLabel <- paste0(round(breaks[uVals], 4), "-",
                        round(breaks[uVals + 1], 4))

  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_categorise", list(breaks = breaks, n = n, layer = out_layer))
  out <- msk_set(obj, out_layer, temp, prov)

  # attach categories
  out@categories[[out_layer]] <- list(gid = uVals, val = groupsLabel)
  out
}
