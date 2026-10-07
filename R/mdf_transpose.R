#' Transpose a mosaik
#'
#' Mirror the grid along its main diagonal, swapping rows and columns. Cell
#' \code{(r, c)} becomes \code{(c, r)}.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param layer [`character(1)`][character]\cr the layer to transpose. If
#'   \code{NULL} (default), all layers are transposed.
#' @return A mosaik with swapped dimensions and transposed cell values.
#' @examples
#' # the map mirrored along its diagonal
#' msk_vis(landscape, .layer("cover"))
#' msk_vis(mdf_transpose(landscape), .layer("cover"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertCharacter
#' @export

mdf_transpose <- function(obj = NULL,
                          layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layer, null.ok = TRUE)

  dims  <- obj@dims   # c(ncols, nrows)
  ext   <- obj@extent
  ncols <- dims[1]
  nrows <- dims[2]

  target_layers <- if(is.null(layer)) names(obj@layers) else layer

  new_layers <- obj@layers
  for(nm in target_layers){
    vals <- msk_pull(obj, nm)
    out <- numeric(length(vals))
    for(r in seq_len(nrows)){
      for(c in seq_len(ncols)){
        old_idx <- (r - 1L) * ncols + c
        new_idx <- (c - 1L) * nrows + r
        out[new_idx] <- vals[old_idx]
      }
    }
    new_layers[[nm]] <- out
  }

  # swap dims and adjust extent
  new_dims <- c(nrows, ncols)
  cx <- (ext[1] + ext[2]) / 2
  cy <- (ext[3] + ext[4]) / 2
  half_w <- (ext[2] - ext[1]) / 2
  half_h <- (ext[4] - ext[3]) / 2
  new_ext <- c(cx - half_h, cx + half_h, cy - half_w, cy + half_w)

  # build output ----

  .update_mosaik(obj,
                 extent = new_ext,
                 dims = new_dims,
                 layers = new_layers,
                 global = list(),
                 step = step)
}
