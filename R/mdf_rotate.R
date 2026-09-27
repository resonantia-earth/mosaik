#' Rotate a mosaik by 90-degree increments
#'
#' Rotate the grid clockwise by 90, 180, or 270 degrees. Extent and dimensions
#' are adjusted accordingly.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param angle [integerish(1)][integer]\cr rotation angle in degrees. Must be
#'   one of \code{90}, \code{180}, or \code{270}.
#' @param layer [character(1)][character]\cr the layer to rotate. If
#'   \code{NULL} (default), all layers are rotated.
#' @return A mosaik with rotated cell values and adjusted extent/dimensions.
#' @examples
#' mdf_rotate(landscape, angle = 90)
#' mdf_rotate(landscape, angle = 180)
#' mdf_rotate(landscape, angle = 270)
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

mdf_rotate <- function(obj = NULL,
                       angle = 90L,
                       layer = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = as.integer(angle), choices = c(90L, 180L, 270L))
  assertCharacter(x = layer, null.ok = TRUE)
  angle <- as.integer(angle)

  dims  <- obj@dims   # c(ncols, nrows)
  ext   <- obj@extent
  ncols <- dims[1]
  nrows <- dims[2]

  target_layers <- if(is.null(layer)) names(obj@layers) else layer

  new_layers <- obj@layers

  for(nm in target_layers){
    vals <- msk_pull(obj, nm)

    if(angle == 90L){
      # clockwise 90: (r,c) -> (c, nrows-1-r)
      out <- numeric(length(vals))
      for(r in seq_len(nrows)){
        for(c in seq_len(ncols)){
          old_idx <- (r - 1L) * ncols + c
          new_r <- c
          new_c <- nrows - r + 1L
          new_idx <- (new_r - 1L) * nrows + new_c
          out[new_idx] <- vals[old_idx]
        }
      }
      new_layers[[nm]] <- out

    } else if(angle == 180L){
      new_layers[[nm]] <- rev(vals)

    } else {
      # clockwise 270 (= counter-clockwise 90): (r,c) -> (ncols-c, r)
      out <- numeric(length(vals))
      for(r in seq_len(nrows)){
        for(c in seq_len(ncols)){
          old_idx <- (r - 1L) * ncols + c
          new_r <- ncols - c + 1L
          new_c <- r
          new_idx <- (new_r - 1L) * nrows + new_c
          out[new_idx] <- vals[old_idx]
        }
      }
      new_layers[[nm]] <- out
    }
  }

  # adjust extent and dims
  if(angle == 180L){
    new_dims <- dims
    new_ext  <- ext
  } else {
    # 90 or 270: swap cols/rows
    new_dims <- c(nrows, ncols)
    cx <- (ext[1] + ext[2]) / 2
    cy <- (ext[3] + ext[4]) / 2
    half_w <- (ext[2] - ext[1]) / 2
    half_h <- (ext[4] - ext[3]) / 2
    new_ext <- c(cx - half_h, cx + half_h, cy - half_w, cy + half_w)
  }

  # build output ----
  prov <- msk_prov("mdf_rotate", list(angle = angle, layer = layer))

  new_mosaik(
    extent     = new_ext,
    dims       = new_dims,
    layers     = new_layers,
    categories = obj@categories,
    patches    = list(),
    global     = list(),
    crs        = obj@crs,
    provenance = c(obj@provenance, list(prov))
  )
}
