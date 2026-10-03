#' Binarise a mosaik
#'
#' Transform a mosaik so that it has the values 0 and 1.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param thresh [`numeric(1)`][numeric]\cr value above which the cell will be
#'   set to 1, below which it will be set to 0.
#' @param match [`numeric(.)`][numeric]\cr one or more values which will be set
#'   to 1, while the remaining values will be set to 0.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}.
#' @family operators to modify cell values
#' @examples
#' # canopy above 10 m, and three classes of the land cover
#' m <- landscape |>
#'   mdf_binarise(thresh = 10, layer = "canopy", add = "above_10") |>
#'   mdf_binarise(match = c(31, 41, 44), layer = "cover", add = "classes")
#' msk_vis(m, .layer("canopy"), .layer("above_10"))
#' msk_vis(m, .layer("cover"), .layer("classes"))
#' @importFrom checkmate assertClass assertNumber assertNumeric assertCharacter
#' @export

mdf_binarise <- function(obj = NULL,
                         thresh = NULL,
                         match = NULL,
                         layer = NULL,
                         add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumber(x = thresh, null.ok = TRUE)
  assertNumeric(x = match, null.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  uVals <- unique(vals)
  if(is.null(thresh)){
    thresh <- max(uVals, na.rm = TRUE)
  } else {
    if(thresh < min(uVals, na.rm = TRUE) || thresh > max(uVals, na.rm = TRUE)){
      stop("please provide a value for 'thresh' within the range of the values of 'obj'.")
    }
  }

  if(!is.null(match)){
    temp <- as.integer(vals %in% match)
  } else {
    temp <- binariseCpp(vals = vals, thresh = thresh)
    temp <- as.integer(as.vector(temp))
  }

  # build output ----
  # binarising changes the layer's kind (-> binary mask); drop any prior role
  .update_mosaik(obj, values = temp, keep = FALSE, step = step)
}
