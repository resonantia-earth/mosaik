#' Determine the centroid of patches
#'
#' The centroid is the average location of all cells of a foreground patch.
#' @param obj [mosaik]\cr the mosaik (with patches) to modify.
#' @param background [integerish(1)][integer]\cr the value any cell with value
#'   NA should have.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, in which the centroid
#'   of each foreground patch has the value of the patch number and all other
#'   cells have the value \code{background}.
#' @examples
#' forest <- mdf_binarise(landscape, match = 41, layer = "cover")
#' forest <- mdf_componentise(forest)
#' mdf_centroid(forest)
#' @family operators to determine objects
#' @importFrom checkmate assertClass assertIntegerish assertCharacter
#' @export

mdf_centroid <- function(obj = NULL,
                         background = NA,
                         layer = NULL,
                         add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertIntegerish(x = background)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  res <- msk_res(obj)

  # build coordinate vectors
  xs <- seq(obj@extent[1] + res[1]/2, by = res[1], length.out = dims[1])
  ys <- seq(obj@extent[3] + res[2]/2, by = res[2], length.out = dims[2])
  coords_x <- rep(xs, times = dims[2])
  coords_y <- rep(ys, each = dims[1])

  # body ----
  dat <- data.frame(val = vals, x = coords_x, y = coords_y)
  datNNA <- dat[!is.na(dat$val),]

  # determine centroids by averaging all cell coordinates per patch
  theMeans <- stats::aggregate(cbind(x, y) ~ val, data = datNNA, FUN = mean)

  # round to nearest cell center
  theMeans$x <- round((theMeans$x - res[1]/2) / res[1]) * res[1] + res[1]/2
  theMeans$y <- round((theMeans$y - res[2]/2) / res[2]) * res[2] + res[2]/2

  # create output grid
  temp <- rep(background, length(vals))
  for(i in seq_len(nrow(theMeans))){
    col_idx <- which.min(abs(xs - theMeans$x[i]))
    row_idx <- which.min(abs(ys - theMeans$y[i]))
    cell_idx <- (row_idx - 1) * dims[1] + col_idx
    temp[cell_idx] <- theMeans$val[i]
  }

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_centroid", list(background = background, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
