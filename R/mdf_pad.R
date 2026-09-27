#' Add or remove cells at the boundary of a mosaik
#'
#' Expand or shrink the grid by adding or removing rows and columns at the
#' edges. Added cells are filled with \code{value}.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param width [integerish(1)][integer]\cr number of cells to add (positive)
#'   or remove (negative) on each side specified by \code{sides}.
#' @param sides [character][character]\cr which sides to pad. Any combination
#'   of \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"}. Defaults
#'   to all four.
#' @param value [numeric(1)][numeric]\cr fill value for added cells. Defaults
#'   to \code{NA}.
#' @param layer [character(1)][character]\cr the layer to modify. If
#'   \code{NULL} (default), all layers are padded.
#' @return A mosaik with modified extent and dimensions.
#' @examples
#' mdf_pad(landscape, width = 5)
#' mdf_pad(landscape, width = 5, sides = c("left", "right"), value = 0)
#' mdf_pad(landscape, width = -3)
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertIntegerish assertSubset assertCharacter
#' @export

mdf_pad <- function(obj = NULL,
                    width = 1L,
                    sides = c("left", "right", "top", "bottom"),
                    value = NA,
                    layer = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertIntegerish(x = width, len = 1, any.missing = FALSE)
  assertSubset(x = sides, choices = c("left", "right", "top", "bottom"))
  assertCharacter(x = layer, null.ok = TRUE)
  width <- as.integer(width)

  dims    <- obj@dims   # c(ncols, nrows)
  ext     <- obj@extent
  res     <- msk_res(obj)
  ncols   <- dims[1]
  nrows   <- dims[2]

  target_layers <- if(is.null(layer)) names(obj@layers) else layer

  if(width > 0){
    # --- add cells ---
    pad_l <- if("left"   %in% sides) width else 0L
    pad_r <- if("right"  %in% sides) width else 0L
    pad_t <- if("top"    %in% sides) width else 0L
    pad_b <- if("bottom" %in% sides) width else 0L

    new_ncols <- ncols + pad_l + pad_r
    new_nrows <- nrows + pad_t + pad_b
    new_ext <- c(ext[1] - pad_l * res[1],
                 ext[2] + pad_r * res[1],
                 ext[3] - pad_b * res[2],
                 ext[4] + pad_t * res[2])

    new_layers <- obj@layers
    for(nm in target_layers){
      vals <- msk_pull(obj, nm)
      out <- rep(value, new_ncols * new_nrows)
      for(r in seq_len(nrows)){
        src_start <- (r - 1L) * ncols + 1L
        dst_start <- (r - 1L + pad_t) * new_ncols + pad_l + 1L
        out[dst_start:(dst_start + ncols - 1L)] <- vals[src_start:(src_start + ncols - 1L)]
      }
      new_layers[[nm]] <- out
    }

  } else if(width < 0){
    # --- remove cells ---
    trim <- abs(width)
    trim_l <- if("left"   %in% sides) trim else 0L
    trim_r <- if("right"  %in% sides) trim else 0L
    trim_t <- if("top"    %in% sides) trim else 0L
    trim_b <- if("bottom" %in% sides) trim else 0L

    new_ncols <- ncols - trim_l - trim_r
    new_nrows <- nrows - trim_t - trim_b
    if(new_ncols < 1L || new_nrows < 1L){
      stop("Trimming removes all cells.")
    }

    new_ext <- c(ext[1] + trim_l * res[1],
                 ext[2] - trim_r * res[1],
                 ext[3] + trim_b * res[2],
                 ext[4] - trim_t * res[2])

    new_layers <- obj@layers
    for(nm in target_layers){
      vals <- msk_pull(obj, nm)
      out <- numeric(new_ncols * new_nrows)
      for(r in seq_len(new_nrows)){
        src_start <- (r - 1L + trim_t) * ncols + trim_l + 1L
        dst_start <- (r - 1L) * new_ncols + 1L
        out[dst_start:(dst_start + new_ncols - 1L)] <- vals[src_start:(src_start + ncols - trim_l - trim_r - 1L)]
      }
      new_layers[[nm]] <- out
    }

  } else {
    return(obj)
  }

  # build output ----
  prov <- msk_prov("mdf_pad", list(width = width, sides = sides,
                                      value = value, layer = layer))

  new_mosaik(
    extent     = new_ext,
    dims       = c(new_ncols, new_nrows),
    layers     = new_layers,
    categories = obj@categories,
    patches    = list(),
    global     = list(),
    crs        = obj@crs,
    provenance = c(obj@provenance, list(prov))
  )
}
