#' Measure distance
#'
#' Measure how far apart the patches of a layer are. Without a cost layer,
#' every cell counts the same, and the distance is measured in cells. With a
#' cost layer, every cell counts with the cost of crossing it, so the distance
#' says how hard it is to get from one patch to another.
#'
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param cost [`character(1)`][character]\cr the name of a layer with the cost
#'   of crossing each cell; \code{NA} marks a cell that cannot be crossed. If
#'   \code{NULL} (default), every cell costs the same.
#' @param routing [`character(1)`][character]\cr which path between two patches
#'   is measured: \code{"straight"} is the direct line, \code{"cheapest"}
#'   (default) the path with the lowest total cost.
#' @param layer [`character(1)`][character]\cr the layer whose patches are
#'   measured. Its patches must have been numbered with
#'   \code{\link{mdf_componentise}} first. Defaults to the first layer.
#' @return The input mosaik with the distances added to the patch results of
#'   \code{layer} (see \code{\link{msk_patches}}).
#' @details The result is called \code{distance}, so that \code{\link{msr}}
#'   reads it as \code{distance.patch}. With a cost layer, it takes the name of
#'   that layer instead: after \code{cost = "friction"}, it is
#'   \code{friction.patch}.
#'
#'   With a cost layer, the distance along a path is the sum of the costs of
#'   all cells on that path. To make paths avoid a class, give that class the
#'   cost \code{NA}, for example with \code{\link{mdf_replace}}. If two patches
#'   cannot be reached from each other, their distance is \code{NA}, not 0.
#'
#'   The result is one matrix for each class. It has one row and one column
#'   for each patch of that class, and each value is the distance between two
#'   of these patches. The distance of a patch to itself is \code{Inf}, so
#'   that \code{min(distance.patch)} in \code{\link{msr}} gives the distance to
#'   the nearest other patch. For a patch that touches the map border, all
#'   distances from that patch are \code{NA}, because another patch beyond the
#'   map could be closer to it. The distances from the other patches to this
#'   patch are still measured.
#'
#'   The distance of every cell to something, such as the edge of its patch,
#'   is a layer; \code{\link{mdf_distance}} computes it.
#' @examples
#' # the forest patches, numbered first
#' f <- mdf_filter(landscape, cover == 47, add = "forest")
#' f <- mdf_componentise(f, connectivity = 8L, layer = "forest", add = "patch")
#'
#' # the distances between the forest patches
#' m <- msr_distance(f, layer = "forest")
#' msk_patches(m, layer = "forest")$distance[["1"]]
#'
#' # the distances over a cost layer built from the classes
#' m <- mdf_replace(f, old = c(21, 24, 47), new = c(5, 2, 1),
#'                  layer = "cover", add = "friction")
#' m <- msr_distance(m, cost = "friction", layer = "forest")
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_distance <- function(obj = NULL, cost = NULL, routing = "cheapest",
                         layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = routing, choices = c("straight", "cheapest"))
  assertCharacter(x = cost, len = 1, null.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)

  if(is.null(layer)) layer <- names(obj@layers)[1]
  if(!(layer %in% names(obj@layers))){
    stop("layer '", layer, "' not found in 'obj'.")
  }
  if(!is.null(cost) && !(cost %in% names(obj@layers))){
    stop("cost layer '", cost, "' not found in 'obj'.")
  }

  # the result is named after what it measures, not after this function
  metric <- if(is.null(cost)) "distance" else cost

  dims <- obj@dims
  vals <- msk_pull(obj, layer)


  surface <- if(is.null(cost)) NULL else msk_pull(obj, cost)

  # the patches numbered by mdf_componentise, one matrix per class
  pt <- .patches_of(obj, layer)
  cc <- pt$ids
  uVals <- sort(unique(pt$class))

  dist_matrices <- list()

  for(i in seq_along(uVals)){

    # this class's patches, in the order of the patch record
    patch_ids <- pt$patch[pt$class == uVals[i]]
    n_patches <- length(patch_ids)

    # the diagonal stays Inf: a patch has no distance to itself, and this
    # keeps min(distance.patch) reading as the nearest neighbour
    if(n_patches <= 1){
      mat <- matrix(Inf, nrow = n_patches, ncol = n_patches)
      if(n_patches == 1){
        rownames(mat) <- colnames(mat) <- as.character(patch_ids)
      }
      dist_matrices[[as.character(uVals[i])]] <- mat
      next
    }

    patch_cells <- lapply(patch_ids, function(pid) which(!is.na(cc) & cc == pid))

    # diagonal stays Inf; unreachable pairs become NA
    mat <- matrix(Inf, nrow = n_patches, ncol = n_patches)
    rownames(mat) <- colnames(mat) <- as.character(patch_ids)

    for(p in seq_len(n_patches)){

      if(is.null(cost) && routing == "straight"){
        # geometric distance needs no path: the transform gives the nearest
        # edge distance directly
        src <- rep(0, length(vals))
        src[patch_cells[[p]]] <- 1
        d <- sqrt(distanceCpp(vals = src, nrow = dims[2], ncol = dims[1],
                              method = "euclidean"))
        for(q in seq_len(n_patches)){
          if(p == q) next
          mat[p, q] <- min(d[patch_cells[[q]]])
        }
        next
      }

      cellCost <- if(is.null(cost)) rep(1, length(vals)) else surface

      if(routing == "straight"){
        for(q in seq_len(n_patches)){
          if(p == q) next
          mat[p, q] <- .straight_path(patch_cells[[p]], patch_cells[[q]],
                                      cellCost, dims)
        }
        next
      }

      # the accumulated cost of the cheapest path, entering the target patch
      # at its cheapest cell
      cd <- costDistanceCpp(cost = cellCost, from = patch_cells[[p]],
                            nrow = dims[2], ncol = dims[1], diagonal = TRUE)
      for(q in seq_len(n_patches)){
        if(p == q) next
        reached <- cd$dist[patch_cells[[q]]]
        reached <- reached[is.finite(reached)]
        mat[p, q] <- if(length(reached) == 0) NA_real_ else min(reached)
      }
    }

    dist_matrices[[as.character(uVals[i])]] <- mat
  }

  # a patch the map border cuts may lie closer to others beyond it; the way
  # from a whole patch to the part of it on the map is still measured
  cut <- as.character(pt$patch[pt$clipped])
  dist_matrices <- lapply(dist_matrices, function(mat){
    mat[rownames(mat) %in% cut, ] <- NA
    mat
  })

  obj@patches[[layer]][[metric]] <- dist_matrices

  obj <- .update_mosaik(obj, step = step)

  return(obj)
}

#' Sum the costs along the direct line between two patches
#'
#' Walks the straight line between the closest pair of cells of two patches and
#' sums the costs of the cells it crosses, whether or not a cheaper detour
#' exists.
#'
#' @param fromCells [`integer(.)`][integer]\cr cell indices of the source patch.
#' @param toCells [`integer(.)`][integer]\cr cell indices of the target patch.
#' @param cellCost [`numeric(.)`][numeric]\cr the cost of each cell.
#' @param dims [`integer(2)`][integer]\cr grid dimensions, columns then rows.
#' @return A single numeric value, \code{NA} when the line crosses a cell that
#'   cannot be crossed.
#' @keywords internal

.straight_path <- function(fromCells, toCells, cellCost, dims){

  ncol <- dims[1]
  fr <- (fromCells - 1) %/% ncol; fc <- (fromCells - 1) %% ncol
  tr <- (toCells - 1) %/% ncol;   tc <- (toCells - 1) %% ncol

  # the closest pair of cells defines the line
  d2 <- outer(fr, tr, function(a, b) (a - b)^2) +
    outer(fc, tc, function(a, b) (a - b)^2)
  best <- which(d2 == min(d2), arr.ind = TRUE)[1, ]

  r0 <- fr[best[1]]; c0 <- fc[best[1]]
  r1 <- tr[best[2]]; c1 <- tc[best[2]]

  # sample the line at every cell it passes through
  n <- max(abs(r1 - r0), abs(c1 - c0)) + 1
  rows <- round(seq(from = r0, to = r1, length.out = n))
  cols <- round(seq(from = c0, to = c1, length.out = n))
  path <- rows * ncol + cols + 1

  vals <- cellCost[path]
  # a cell that cannot be crossed makes the direct route unavailable
  if(any(is.na(vals))) return(NA_real_)

  sum(vals)
}
