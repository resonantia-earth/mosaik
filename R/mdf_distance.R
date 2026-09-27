#' Calculate the distance map for a mosaik
#'
#' Replace each cell value with its distance to the nearest source, where the
#' source is one edge of a binary layer, a set of labelled target points, or
#' explicit coordinates.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param source [character(1)][character]\cr what to measure distance from:
#'   \code{"foreground"} (default) or \code{"background"} for the two edges of a
#'   binary layer, or the name of a layer holding target points with matching
#'   patch IDs. Ignored when \code{coords} is provided. See Details.
#' @param coords [matrix][matrix]\cr optional two-column matrix of source
#'   coordinates (x, y) in CRS space. When provided, \code{source} is ignored
#'   and \code{obj} does not need to be binary.
#' @param snap [logical(1)][logical]\cr only used when \code{coords} is
#'   provided. If \code{TRUE} (default), source coordinates are snapped to the
#'   nearest cell centroid and the fast Meijster algorithm is used. If
#'   \code{FALSE}, exact Euclidean distances from the original coordinates are
#'   computed (slower but precise).
#' @param method [character(1)][character]\cr the distance measure to
#'   calculate. Either \code{"euclidean"} (default), \code{"manhattan"} or
#'   \code{"chessboard"}. Only used with the Meijster algorithm (binary layer
#'   or \code{snap = TRUE}).
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, where cell values have
#'   been replaced with the distance to the nearest source. Distances are in
#'   cell units when using the Meijster algorithm, or in CRS units when using
#'   exact coordinates (\code{snap = FALSE}).
#' @details
#'   The three modes of \code{source}:
#'   \describe{
#'     \item{\code{"foreground"}}{distance from each background cell to the
#'       nearest foreground cell. Input must be binary.}
#'     \item{\code{"background"}}{distance from each foreground cell to the
#'       nearest background cell (internal edge distance). Input must be binary.}
#'     \item{a layer name}{distance from each cell to the target point sharing
#'       its value. \code{layer} holds patch IDs (e.g. from
#'       \code{\link{mdf_componentise}}) and the named \code{source} layer holds
#'       sparse points with matching IDs (e.g. from \code{\link{mdf_centroid}});
#'       cells with no matching target get \code{NA}.}
#'   }
#' @examples
#' # distance from foreground cells in a binary layer
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' mdf_distance(forest)
#' mdf_distance(forest, method = "manhattan")
#'
#' # distance from each foreground cell to nearest background (edge distance)
#' mdf_distance(forest, source = "background")
#'
#' # distance from each cell to its own patch centroid (GYRATE workflow)
#' forest |>
#'   mdf_componentise(add = "patches") |>
#'   mdf_centroid(layer = "patches", add = "centroids") |>
#'   mdf_distance(source = "centroids", layer = "patches", add = "dist")
#'
#' # distance from explicit coordinates
#' mdf_distance(landscape, coords = matrix(c(30, 28), ncol = 2),
#'              layer = "cover")
#' @references Meijster, A., Roerdink, J.B.T.M., Hesselink, W.H., 2000. A
#'   general algorithm for computing distance transforms in linear time, in:
#'   Goutsias, J., Vincent, L., Bloomberg, D.S. (Eds.), Mathematical Morphology
#'   and Its Applications to Image and Signal Processing. Springer, pp. 331-340.
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertString assertCharacter assertLogical assertMatrix
#' @export

mdf_distance <- function(obj = NULL,
                         source = "foreground",
                         coords = NULL,
                         snap = TRUE,
                         method = "euclidean",
                         layer = NULL,
                         add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertString(x = source)
  assertChoice(x = method, choices = c("euclidean", "manhattan", "chessboard"))
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)
  assertLogical(x = snap, len = 1, any.missing = FALSE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  dims <- obj@dims
  ext  <- obj@extent
  res  <- msk_res(obj)

  if(!is.null(coords)){

    assertMatrix(x = coords, ncols = 2, mode = "numeric", min.rows = 1)

    if(snap){
      # snap coords to nearest cell centroids, build binary grid, run Meijster
      binary <- rep(0, prod(dims))
      xs <- seq(ext[1] + res[1] / 2, by = res[1], length.out = dims[1])
      ys <- seq(ext[4] - res[2] / 2, by = -res[2], length.out = dims[2])

      for(i in seq_len(nrow(coords))){
        col_idx <- which.min(abs(xs - coords[i, 1]))
        row_idx <- which.min(abs(ys - coords[i, 2]))
        cell_idx <- (row_idx - 1L) * dims[1] + col_idx
        binary[cell_idx] <- 1
      }

      temp <- distanceCpp(vals = binary, nrow = dims[2], ncol = dims[1],
                          method = method)
      if(method == "euclidean") temp <- sqrt(temp)

    } else {
      # exact distance from coordinates to cell centroids
      if(method != "euclidean"){
        warning("method '", method, "' ignored when snap = FALSE; ",
                "exact mode always uses Euclidean distance.")
      }
      temp <- distanceFromPointsCpp(xmin = ext[1], ymax = ext[4],
                                     res_x = res[1], res_y = res[2],
                                     ncol = dims[1], nrow = dims[2],
                                     coords = coords)
    }

  } else if(source %in% c("foreground", "background")){
    # binary distance transform
    vals <- msk_pull(obj, layer)
    if(!isBinaryCpp(vals = vals)){
      stop("'obj' is not binary, please run 'mdf_binarise()' first or provide 'coords'.")
    }

    if(source == "background"){
      vals <- as.integer(vals == 0)
    }

    temp <- distanceCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                        method = method)
    if(method == "euclidean") temp <- sqrt(temp)

  } else {
    # layer-based mode: source is a layer name with target points
    if(!(source %in% names(obj@layers))){
      stop("source layer '", source, "' not found in the mosaik.")
    }

    vals <- msk_pull(obj, layer)
    targets <- msk_pull(obj, source)
    ncols <- dims[1]
    nrows <- dims[2]

    # cell coordinates in grid units (0.5-indexed)
    col_coords <- rep(seq_len(ncols), times = nrows) - 0.5
    row_coords <- rep(seq_len(nrows), each = ncols) - 0.5

    # build lookup: patch ID -> target cell index
    target_cells <- which(!is.na(targets) & targets != 0)
    if(length(target_cells) == 0){
      stop("source layer '", source, "' has no non-NA/non-zero target cells.")
    }
    target_ids <- targets[target_cells]
    target_x <- col_coords[target_cells]
    target_y <- row_coords[target_cells]

    # for each target ID, store its coordinates
    tgt_lookup <- split(seq_along(target_cells),  target_ids)

    temp <- rep(NA_real_, length(vals))

    for(id_str in names(tgt_lookup)){
      id <- as.numeric(id_str)
      member_cells <- which(vals == id)
      if(length(member_cells) == 0) next

      # target point(s) for this ID
      tidx <- tgt_lookup[[id_str]]
      tx <- target_x[tidx]
      ty <- target_y[tidx]

      for(ci in member_cells){
        # distance to nearest target with this ID
        dists <- sqrt((col_coords[ci] - tx)^2 + (row_coords[ci] - ty)^2)
        temp[ci] <- min(dists)
      }
    }
  }

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_distance", list(source = source, coords = coords,
                                           snap = snap, method = method,
                                           layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
