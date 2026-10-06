#' Fill enclosed holes in foreground patches
#'
#' Turn background cells that are wholly surrounded by foreground into
#' foreground, closing the holes of a patch while leaving its outer shape
#' unchanged. Background that reaches the edge of the grid is not a hole and
#' stays background.
#'
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param connectivity [`integerish(1)`][integer]\cr how background cells
#'   connect when deciding what is enclosed: \code{4} (default, rook) or
#'   \code{8} (queen). Under \code{4}, a hole that escapes to the edge only
#'   through a diagonal counts as enclosed and is filled.
#' @param layer [`character(1)`][character]\cr the binary layer in \code{obj}
#'   to fill. Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik in which the holes of \code{layer} are 1.
#' @details
#'   The layer must be binary: 1 is foreground, 0 and \code{NA} are background.
#'   To fill the holes of one class of a categorical map, binarise it first
#'   with \code{\link{mdf_binarise}}.
#'
#'   The background is split into connected areas, and an area is a hole when
#'   none of its cells lie on the edge of the grid. Because the test is
#'   connectivity and not a moving window, a hole is filled whatever its size.
#'   A morphological close (\code{\link{mdf_dilate}} then
#'   \code{\link{mdf_erode}}) only closes gaps up to the size of its
#'   structuring element, and \code{\link{mdf_dilate}} alone grows the patch
#'   outward as well.
#' @family operators to morphologically modify a raster
#' @examples
#' # the forest with its holes filled, and the holes alone
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_fill(layer = "forest", add = "filled") |>
#'   mdf_filter(filled == 1 & forest == 0, add = "holes")
#' msk_vis(m, .layer("forest"), .layer("filled"), .layer("holes"))
#' @importFrom checkmate assertClass assertCharacter assertChoice
#' @export

mdf_fill <- function(obj = NULL,
                     connectivity = 4L,
                     layer = NULL,
                     add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = connectivity, choices = c(4L, 8L))
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  if(!isBinaryCpp(vals = vals)){
    stop("layer '", layer, "' is not binary; make it binary with 'mdf_filter()' first.",
         call. = FALSE)
  }
  dims <- obj@dims
  ncols <- dims[1]
  nrows <- dims[2]

  # body ----
  # background mask (1 = background, NA = foreground) so componentsCpp labels
  # only the background regions; foreground is a single NA sea between them.
  is_bg <- is.na(vals) | vals == 0
  bg <- rep(NA_real_, length(vals))
  bg[is_bg] <- 1

  cc <- componentsCpp(vals = bg, nrow = nrows, ncol = ncols,
                      connectivity = connectivity)

  # a background component is enclosed iff none of its cells sit on the grid edge
  if(nrows == 1 || ncols == 1){
    edge <- seq_len(nrows * ncols)
  } else {
    edge <- unique(c(seq_len(ncols),                          # top row
                     (nrows - 1) * ncols + seq_len(ncols),    # bottom row
                     (seq_len(nrows) - 1) * ncols + 1,        # left column
                     seq_len(nrows) * ncols))                 # right column
  }
  edge_ids <- unique(cc[edge])
  edge_ids <- edge_ids[!is.na(edge_ids)]
  enclosed <- !is.na(cc) & !(cc %in% edge_ids)

  temp <- vals
  temp[enclosed] <- 1

  # build output ----
  .update_mosaik(obj, values = temp, step = step)
}
