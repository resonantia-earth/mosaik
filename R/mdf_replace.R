#' Replace values in a mosaik
#'
#' Replace a set of values with new values.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param old [numeric(.)][numeric]\cr values to be substituted.
#' @param new [numeric(.)][numeric]\cr values to substitute with.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}.
#' @examples
#' mdf_replace(landscape, old = 47, new = 99, layer = "cover")
#' mdf_replace(landscape, old = c(21, 24), new = c(100, 200),
#'             layer = "cover")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertNumeric assertCharacter
#' @export

mdf_replace <- function(obj = NULL,
                        old = NULL,
                        new = NULL,
                        layer = NULL,
                        add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  # NA in 'old' is allowed: it addresses the cells that carry no value, which is
  # how a masked layer (mdf_filter) is turned back into a 0/1 mask. The body
  # branches on all(is.na(old)) for exactly this case.
  assertNumeric(x = old, min.len = 1)
  # NA in 'new' is allowed: replacing a value with NA marks those cells as
  # excluded, e.g. an impassable class for msr_cost
  assertNumeric(x = new, min.len = 1)
  if(length(old) != length(new)){
    newValues <- rep(new, length.out = length(old))
  } else{
    newValues <- new
  }
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  if(all(is.na(old))){
    temp <- vals
    temp[is.na(temp)] <- new
  } else {
    temp <- vals
    for(i in seq_along(old)){
      temp[temp == old[i]] <- newValues[i]
    }
  }

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_replace", list(old = old, new = new, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
