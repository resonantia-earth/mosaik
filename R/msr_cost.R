#' Measure the cost between patches
#'
#' Measure what it costs to get from one patch to another. Distance is the
#' member of that family measured in metres; any other cost is expressed by
#' handing in a surface of per-cell traversal costs.
#'
#' @param obj [mosaik]\cr the mosaik to measure.
#' @param scale [character(1)][character]\cr \code{"patch"} (default) for
#'   pairwise costs between patches, or \code{"cell"} for a per-cell cost
#'   surface stored as a layer.
#' @param cost [character(1)][character]\cr the name of a layer holding
#'   per-cell traversal costs. If \code{NULL} (default), cost is geometric
#'   distance in metres.
#' @param routing [character(1)][character]\cr how the path between two patches
#'   is chosen: \code{"straight"} follows the direct line, \code{"cheapest"}
#'   (default) searches for the path of least accumulated cost. Only used at
#'   \code{scale = "patch"}.
#' @param accumulate [character(1)][character]\cr how the cell values along the
#'   path are combined into one number: \code{"sum"} (default), \code{"max"},
#'   \code{"min"}, \code{"product"} or \code{"mean"}. Only used at
#'   \code{scale = "patch"}.
#' @param layer [character(1)][character]\cr the layer holding the patches to
#'   measure between. Defaults to the first layer.
#' @return The input mosaik with cost information added. At
#'   \code{scale = "patch"} a named list of matrices, one per class, is written
#'   to \code{@patches}; at \code{scale = "cell"} an internal layer is added.
#'   Both are named after what the surface measures rather than after this
#'   function: \code{distance} when \code{cost} is \code{NULL}, otherwise the
#'   name of the \code{cost} layer. So a friction surface yields
#'   \code{friction.patch} and \code{friction.cell} in \code{\link{msr}()}.
#' @details \code{NA} in the cost surface means impassable, following the
#'   convention of the other \code{msr_*} functions. Avoiding a class is
#'   therefore done by setting it to \code{NA} in the surface, e.g. with
#'   \code{\link{mdf_replace}()}; there is no separate argument for it. A pair
#'   of patches with no route between them is \code{NA}, which is distinct from
#'   a cost of zero.
#'
#'   \code{accumulate} is applied to the cells along the chosen path.
#'   \code{"sum"} suits a surface holding a per-cell price, while \code{"max"}
#'   is what an already-accumulated surface needs: the values of a distance map
#'   grow away from their source, so the largest one along the path is the value
#'   at its end. \code{"product"} suits per-cell passage probabilities.
#'
#'   Note that \code{routing} and \code{accumulate} are independent.
#'   \code{"cheapest"} finds the path minimising the accumulated cost, so
#'   combining it with \code{accumulate = "max"} reports the worst cell on the
#'   cheapest path — not the path whose worst cell is lowest.
#'
#'   At \code{scale = "cell"} the surface is written under a reserved name
#'   prefixed with an underscore, which is how \code{\link{msr}()} finds it. The
#'   two scales aggregate different things, both to one value per patch:
#'   \code{max(distance.patch)} is the cost to the most distant other patch,
#'   whereas \code{max(distance.cell)} is the most expensive cell within a
#'   patch.
#'
#'   \code{scale = "cell"} expects \code{layer} to be \strong{binary}: the
#'   surface measures the cost of reaching each cell from the background, so it
#'   needs a foreground to measure into. Isolate the class of interest with
#'   \code{\link{mdf_binarise}()} first. Reducing that surface per patch with
#'   \code{\link{msr}()} then needs a layer of patch IDs to group by — see the
#'   examples — because \code{msr()} groups cell values by the values of the
#'   layer it is given, and a binary layer has only two of those.
#' @examples
#' # pairwise distance between patches, in metres
#' m <- msr_cost(landscape)
#' m@patches$distance[["47"]]
#'
#' # cost over a friction surface built from the classes
#' m <- mdf_replace(landscape, old = c(21, 24, 47), new = c(5, 2, 1),
#'                  add = "friction")
#' m <- msr_cost(m, cost = "friction", layer = "cover")
#'
#' # the worst barrier on the cheapest path, rather than the total
#' m <- msr_cost(m, cost = "friction", accumulate = "max", layer = "cover")
#'
#' # a per-cell surface, reduced per patch by msr(). The surface measures how
#' # deep into a patch each cell lies, so it needs one class on its own: a
#' # binary layer to measure into, and its patch IDs to group the result by.
#' f <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest")
#' f <- mdf_componentise(f, layer = "forest", add = "patches")
#' f <- msr_cost(f, scale = "cell", layer = "forest")
#'
#' # the deepest cell of each patch — its distance to the nearest edge
#' f <- msr(f, equation = "max(distance.cell)", label = "core",
#'          layer = "patches")
#' f@patches$core
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_cost <- function(obj, scale = "patch", cost = NULL, routing = "cheapest",
                     accumulate = "sum", layer = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = scale, choices = c("patch", "cell"))
  assertChoice(x = routing, choices = c("straight", "cheapest"))
  assertChoice(x = accumulate, choices = c("sum", "max", "min", "product", "mean"))
  assertCharacter(x = cost, len = 1, null.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)

  if(is.null(layer)) layer <- names(obj@layers)[1]
  if(!(layer %in% names(obj@layers))){
    stop("layer '", layer, "' not found in 'obj'.")
  }
  if(!is.null(cost) && !(cost %in% names(obj@layers))){
    stop("cost layer '", cost, "' not found in 'obj'.")
  }

  # the surface is named after what it measures, not after this function
  metric <- if(is.null(cost)) "distance" else cost

  dims <- obj@dims
  vals <- msk_pull(obj, layer)

  # --- scale = "cell": register the surface for msr() ----
  if(scale == "cell"){

    if(is.null(cost)){
      obj <- mdf_distance(obj, source = "background", layer = layer,
                          add = paste0("_", metric))
    } else {
      surface <- msk_pull(obj, cost)
      src <- which(!is.na(vals) & vals != 0)
      cd <- costDistanceCpp(cost = surface, from = src, nrow = dims[2],
                            ncol = dims[1], diagonal = TRUE)
      out <- cd$dist
      out[is.infinite(out)] <- NA_real_
      obj@layers[[paste0("_", metric)]] <- out
    }

    prov <- msk_prov("msr_cost", list(scale = scale, cost = cost,
                                        layer = layer))
    obj@provenance <- c(obj@provenance, list(prov))
    return(obj)
  }

  # --- scale = "patch": pairwise cost between patches ----
  surface <- if(is.null(cost)) NULL else msk_pull(obj, cost)
  uVals <- sort(unique(vals[!is.na(vals)]))

  # one global componentisation, so patch IDs are unique across all classes and
  # index the same node set as msr_adjacency(scale = "patch")
  cc <- componentsCpp(vals = vals, nrow = dims[2], ncol = dims[1])

  cost_matrices <- list()

  for(i in seq_along(uVals)){

    # this class's patches, keeping their global IDs
    class_cc <- cc
    class_cc[vals != uVals[i]] <- NA
    patch_ids <- sort(unique(class_cc[!is.na(class_cc)]))
    n_patches <- length(patch_ids)

    # the diagonal stays Inf: a patch has no cost to itself, and this keeps
    # min(cost.patch) reading as the nearest neighbour
    if(n_patches <= 1){
      mat <- matrix(Inf, nrow = n_patches, ncol = n_patches)
      if(n_patches == 1){
        rownames(mat) <- colnames(mat) <- as.character(patch_ids)
      }
      cost_matrices[[as.character(uVals[i])]] <- mat
      next
    }

    patch_cells <- lapply(patch_ids, function(pid) which(cc == pid))

    # diagonal stays Inf (no cost to self); unreachable pairs become NA
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
                                      cellCost, accumulate, dims)
        }
        next
      }

      cd <- costDistanceCpp(cost = cellCost, from = patch_cells[[p]],
                            nrow = dims[2], ncol = dims[1], diagonal = TRUE)

      for(q in seq_len(n_patches)){
        if(p == q) next
        mat[p, q] <- .accumulate_path(cd, patch_cells[[q]], cellCost,
                                      accumulate)
      }
    }

    cost_matrices[[as.character(uVals[i])]] <- mat
  }

  obj@patches[[metric]] <- cost_matrices

  prov <- msk_prov("msr_cost", list(scale = scale, cost = cost,
                                      routing = routing,
                                      accumulate = accumulate, layer = layer))
  obj@provenance <- c(obj@provenance, list(prov))

  return(obj)
}

#' Combine the cell values along a path
#'
#' Traces the cheapest path back through the predecessor map and reduces the
#' cells it crosses. \code{"sum"} is read off the accumulated surface directly,
#' since that is what Dijkstra already computed.
#'
#' @param cd [list][list]\cr the return value of \code{costDistanceCpp}.
#' @param targetCells [integer(.)][integer]\cr cell indices of the target patch.
#' @param cellCost [numeric(.)][numeric]\cr the per-cell cost surface.
#' @param accumulate [character(1)][character]\cr the reduction to apply.
#' @return A single numeric value, \code{NA} when the target is unreachable.
#' @keywords internal

.accumulate_path <- function(cd, targetCells, cellCost, accumulate){

  reached <- targetCells[is.finite(cd$dist[targetCells])]
  if(length(reached) == 0) return(NA_real_)

  # enter the target patch at its cheapest cell
  endCell <- reached[which.min(cd$dist[reached])]

  # the accumulated surface already holds the summed cost
  if(accumulate == "sum") return(cd$dist[endCell])

  # otherwise the cells actually crossed are needed
  path <- endCell
  cur <- endCell
  while(cd$pred[cur] > 0){
    cur <- cd$pred[cur]
    path <- c(path, cur)
  }

  .reduce_path(cellCost[path], accumulate)
}

#' Combine the cell values along the direct line between two patches
#'
#' Walks the straight line between the closest pair of cells of two patches and
#' reduces the cells it crosses, whether or not a cheaper detour exists.
#'
#' @param fromCells [integer(.)][integer]\cr cell indices of the source patch.
#' @param toCells [integer(.)][integer]\cr cell indices of the target patch.
#' @param cellCost [numeric(.)][numeric]\cr the per-cell cost surface.
#' @param accumulate [character(1)][character]\cr the reduction to apply.
#' @param dims [integer(2)][integer]\cr grid dimensions, columns then rows.
#' @return A single numeric value, \code{NA} when the line crosses an
#'   impassable cell.
#' @keywords internal

.straight_path <- function(fromCells, toCells, cellCost, accumulate, dims){

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
  # an impassable cell on the line makes the direct route unavailable
  if(any(is.na(vals))) return(NA_real_)

  .reduce_path(vals, accumulate)
}

#' Reduce the cell values of a path to a single number
#'
#' @param vals [numeric(.)][numeric]\cr the cell values along the path.
#' @param accumulate [character(1)][character]\cr the reduction to apply.
#' @return A single numeric value, \code{NA} when nothing remains to reduce.
#' @keywords internal

.reduce_path <- function(vals, accumulate){

  vals <- vals[!is.na(vals)]
  if(length(vals) == 0) return(NA_real_)

  switch(accumulate,
         sum = sum(vals),
         max = max(vals),
         min = min(vals),
         mean = mean(vals),
         product = prod(vals))
}
