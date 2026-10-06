#' Replace values in a mosaik
#'
#' Replace a set of values with new values.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param old [`numeric(.)`][numeric]\cr values to be substituted.
#' @param new [`numeric(.)`][numeric]\cr values to substitute with.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}.
#' @examples
#' # class 47 renumbered to 99, and classes 21 and 24 to 100 and 200
#' m <- landscape |>
#'   mdf_replace(old = 47, new = 99, layer = "cover", add = "one") |>
#'   mdf_replace(old = c(21, 24), new = c(100, 200), layer = "cover",
#'               add = "two")
#' msk_vis(m, .layer("cover"), .layer("one"), .layer("two"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertNumeric assertCharacter
#' @export

mdf_replace <- function(obj = NULL,
                        old = NULL,
                        new = NULL,
                        layer = NULL,
                        add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  # NA in 'old' is allowed: it addresses the cells that carry no value, which is
  # how a masked layer (mdf_filter) is turned back into a 0/1 mask. The body
  # branches on all(is.na(old)) for exactly this case.
  assertNumeric(x = old, min.len = 1)
  # NA in 'new' is allowed: replacing a value with NA marks those cells as
  # excluded, e.g. an impassable class for msr_distance
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
    # all values are replaced at once, so a new value that equals a later old
    # value is not replaced a second time
    temp <- vals
    hit <- match(vals, old)
    temp[!is.na(hit)] <- newValues[hit[!is.na(hit)]]
  }

  # build output ----
  .update_mosaik(obj, values = temp, step = step)
}
