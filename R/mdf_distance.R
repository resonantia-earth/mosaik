#' Calculate the distance map for a mosaik
#'
#' Replace each cell value with its distance to the nearest source, where the
#' source is one side of a binary layer or a set of labelled target cells.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param source [`character(1)`][character]\cr what to measure distance from:
#'   \code{"foreground"} (default) or \code{"background"} for the two sides of a
#'   binary layer, or the name of a layer holding target cells with matching
#'   patch numbers. See Details.
#' @param method [`character(1)`][character]\cr the distance measure to
#'   calculate. Either \code{"euclidean"} (default), \code{"manhattan"} or
#'   \code{"chessboard"}. Applies to \code{"foreground"} and
#'   \code{"background"}; distances to target cells are always Euclidean.
#' @param cost [`character(1)`][character]\cr the name of a layer with the cost
#'   of crossing each cell; \code{NA} marks a cell that cannot be crossed. With
#'   a cost layer, each cell gets the lowest total cost of reaching it from
#'   the source, and \code{method} is not used. Only for \code{"foreground"}
#'   and \code{"background"}. If \code{NULL} (default), every cell costs the
#'   same.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik in which each cell carries its distance to the nearest
#'   source, in cells, or with \code{cost} the lowest total cost of reaching
#'   it.
#' @details
#'   The three modes of \code{source}:
#'   \describe{
#'     \item{\code{"foreground"}}{distance from each background cell to the
#'       nearest foreground cell. Input must be binary.}
#'     \item{\code{"background"}}{distance from each foreground cell to the
#'       nearest background cell (internal edge distance). Input must be binary.}
#'     \item{a layer name}{distance from each cell to the target cell sharing
#'       its value. \code{layer} holds patch numbers (e.g. from
#'       \code{\link{mdf_componentise}}) and the named \code{source} layer holds
#'       single cells with matching numbers (e.g. from
#'       \code{\link{mdf_centroid}}); cells with no matching target get
#'       \code{NA}.}
#'   }
#'
#'   To measure the distance from points, lines or polygons, rasterise them
#'   with \code{\link{msk_rasterise}} first and binarise the result.
#' @examples
#' # distance from the forest, outward and inward
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_distance(layer = "forest", add = "outward") |>
#'   mdf_distance(source = "background", layer = "forest", add = "inward")
#' msk_vis(m, .layer("forest"), .layer("outward"), .layer("inward"))
#'
#' # euclidean, manhattan and chessboard distance from the forest
#' m <- m |>
#'   mdf_distance(method = "manhattan", layer = "forest", add = "manhattan") |>
#'   mdf_distance(method = "chessboard", layer = "forest", add = "chessboard")
#' msk_vis(m, .layer("outward"), .layer("manhattan"), .layer("chessboard"))
#'
#' # the cost of reaching each cell from the forest, over a cost layer built
#' # from the classes
#' m <- m |>
#'   mdf_replace(old = c(21, 24, 47), new = c(5, 2, 1), layer = "cover",
#'               add = "friction") |>
#'   mdf_distance(cost = "friction", layer = "forest", add = "reach")
#' msk_vis(m, .layer("friction"), .layer("reach"))
#'
#' # distance from each cell to the centroid of its own patch
#' m <- m |>
#'   mdf_componentise(layer = "forest", add = "patches") |>
#'   mdf_centroid(layer = "patches", add = "centroids") |>
#'   mdf_distance(source = "centroids", layer = "patches", add = "to_centroid")
#' msk_vis(m, .layer("patches"), .layer("to_centroid"))
#' @references Meijster, A., Roerdink, J.B.T.M., Hesselink, W.H., 2000. A
#'   general algorithm for computing distance transforms in linear time, in:
#'   Goutsias, J., Vincent, L., Bloomberg, D.S. (Eds.), Mathematical Morphology
#'   and Its Applications to Image and Signal Processing. Springer, pp. 331-340.
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertString assertCharacter assertChoice
#' @export

mdf_distance <- function(obj = NULL,
                         source = "foreground",
                         method = "euclidean",
                         cost = NULL,
                         layer = NULL,
                         add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertString(x = source)
  assertChoice(x = method, choices = c("euclidean", "manhattan", "chessboard"))
  assertCharacter(x = cost, len = 1, null.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  dims <- obj@dims
  if(!is.null(cost) && !(cost %in% names(obj@layers))){
    stop("cost layer '", cost, "' not found in 'obj'.")
  }

  if(source %in% c("foreground", "background")){
    # binary distance transform
    vals <- msk_pull(obj, layer)
    if(!isBinaryCpp(vals = vals)){
      stop("'obj' is not binary, make it binary with 'mdf_filter()' first.")
    }

    if(source == "background"){
      vals <- as.integer(vals == 0)
    }

    if(is.null(cost)){
      temp <- distanceCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                          method = method)
      if(method == "euclidean") temp <- sqrt(temp)
    } else {
      # the cheapest accumulated cost from any source cell
      cd <- costDistanceCpp(cost = msk_pull(obj, cost),
                            from = which(!is.na(vals) & vals == 1),
                            nrow = dims[2], ncol = dims[1], diagonal = TRUE)
      temp <- cd$dist
      temp[is.infinite(temp)] <- NA_real_
    }

  } else {
    if(!is.null(cost)){
      stop("'cost' works only with source = \"foreground\" or \"background\".")
    }
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
  .update_mosaik(obj, values = temp, step = step)
}
