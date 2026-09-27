#' Fill enclosed holes in foreground patches
#'
#' Turn background cells that are wholly surrounded by foreground into
#' foreground, closing the interior holes of a patch while leaving its outer
#' shape unchanged. A gap counts as a hole only if it is sealed off from the
#' edge of the grid; any background that still reaches the edge stays background.
#'
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param background [numeric(1)][numeric]\cr the value that counts as
#'   background (holes are made of these). \code{NA} cells are always treated as
#'   background too. Default \code{0}.
#' @param value [numeric(1)][numeric]\cr the value written into filled holes. If
#'   \code{NULL} (default), the most common foreground value is used, so a
#'   single-class mask fills with its own class.
#' @param connectivity [integerish(1)][integer]\cr background connectivity for
#'   deciding what is enclosed: \code{4} (default, rook) or \code{8} (queen). Use
#'   \code{4} to treat diagonally-pinched holes as enclosed.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}.
#' @details
#'   The layer is read as binary: cells equal to \code{background} (and any
#'   \code{NA}) are background, everything else is foreground. The background is
#'   split into connected components with \code{componentsCpp}; a component is
#'   \emph{enclosed} exactly when none of its cells lie on the grid's outer edge,
#'   since a component touching the edge has an open path out of the grid. Every
#'   enclosed component is rewritten to \code{value}; edge-connected background
#'   and the outer boundary of each patch are untouched.
#'
#'   Because the test is connectivity, not a moving window, a hole is filled
#'   whatever its size. A morphological close (\code{\link{mdf_dilate}} then
#'   \code{\link{mdf_erode}}) only closes gaps up to the size of its structuring
#'   element, and \code{\link{mdf_dilate}} alone grows the patch outward as well;
#'   \code{mdf_fill} does neither. \code{connectivity} decides borderline holes
#'   that escape only through a diagonal: under rook (4) the diagonal does not
#'   connect, so the hole is sealed and fills; under queen (8) it leaks to the
#'   edge and is kept.
#' @family operators to morphologically modify a raster
#' @examples
#' # fill the interior gaps of forest patches, leaving patch outlines intact
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' mdf_fill(forest)
#' @importFrom checkmate assertClass assertNumber assertCharacter assertChoice
#' @export

mdf_fill <- function(obj = NULL,
                     background = 0,
                     value = NULL,
                     connectivity = 4L,
                     layer = NULL,
                     add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumber(x = background)
  assertNumber(x = value, null.ok = TRUE)
  assertChoice(x = connectivity, choices = c(4L, 8L))
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  ncols <- dims[1]
  nrows <- dims[2]

  # body ----
  # background mask (1 = background, NA = foreground) so componentsCpp labels
  # only the background regions; foreground is a single NA sea between them.
  is_bg <- is.na(vals) | vals == background
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

  # value written into filled holes: explicit, else the modal foreground value
  if(is.null(value)){
    fg <- vals[!is_bg]
    if(length(fg) == 0){
      stop("'obj' has no foreground to fill with; set 'value' explicitly.")
    }
    value <- as.numeric(names(which.max(table(fg))))
  }

  temp <- vals
  temp[enclosed] <- value

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_fill", list(background = background, value = value,
                     connectivity = connectivity, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
