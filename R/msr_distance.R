#' Measure distance
#'
#' Measure how far apart the classes of a layer are. Without a cost layer,
#' every cell counts the same, and the distance is measured in cells or in map
#' units. With a cost layer, every cell counts with the cost of crossing it, so
#' the distance says how hard it is to get from one class to another.
#'
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param unit [`character(1)`][character]\cr \code{"map"} (default, in map
#'   units) or \code{"cells"} (the distance in cells). With a cost layer,
#'   \code{"map"} multiplies the cost of each cell by its width, so a cost per
#'   map unit gives the cost of the whole path.
#' @param cost [`character(1)`][character]\cr the name of a layer with the cost
#'   of crossing each cell; \code{NA} marks a cell that cannot be crossed. If
#'   \code{NULL} (default), every cell costs the same.
#' @param routing [`character(1)`][character]\cr which path between two classes
#'   is measured: \code{"straight"} is the direct line, \code{"cheapest"}
#'   (default) the path with the lowest total cost.
#' @param name [`character(1)`][character]\cr the name the distances are stored
#'   under, by default \code{"distance"}. Give distances over a cost layer a
#'   name of their own, so that both can be stored. The name must not be the
#'   name of a layer and must not contain \code{.} or \code{_}, because
#'   \code{\link{msr}} reads these as the start of the focus and of the layer.
#' @param layer [`character(1)`][character]\cr the layer whose classes are
#'   measured. Defaults to the first layer.
#' @return The input mosaik with the distances added to the class table of
#'   \code{layer} (see \code{\link{msk_table}}).
#' @details The result is a class by class matrix, in the order of the class
#'   table, and each value is the distance between the nearest cells of two
#'   classes. On a layer of patch numbers from \code{\link{mdf_componentise}},
#'   these are the distances between patches. The distance of a class to
#'   itself is \code{Inf}, so that \code{min(distance.others)} in
#'   \code{\link{msr}} gives the distance to the nearest other class.
#'
#'   With a cost layer, the distance along a path is the sum of the costs of
#'   all cells on that path. To make paths avoid a class, give that class the
#'   cost \code{NA}, for example with \code{\link{mdf_replace}}. If two classes
#'   cannot be reached from each other, their distance is \code{NA}, not 0.
#'
#'   The distance of every cell to something, such as the edge of its patch,
#'   is a layer; \code{\link{mdf_distance}} computes it.
#' @examples
#' # the forest patches, numbered first
#' f <- mdf_filter(landscape, cover == 47, add = "forest") |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch")
#'
#' # the distances between the forest patches
#' m <- msr_distance(f, layer = "patch")
#' msk_table(m, layer = "patch")$distance
#'
#' # the distances over a cost layer built from the classes
#' m <- mdf_replace(m, old = c(21, 24, 47), new = c(5, 2, 1),
#'                  layer = "cover", add = "friction") |>
#'   msr_distance(cost = "friction", name = "effort", layer = "patch")
#' msk_table(m, layer = "patch")$effort
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_distance <- function(obj = NULL, cost = NULL, routing = "cheapest",
                         unit = "map", name = "distance", layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = routing, choices = c("straight", "cheapest"))
  assertChoice(x = unit, choices = c("cells", "map"))
  assertCharacter(x = cost, len = 1, null.ok = TRUE)
  assertCharacter(x = name, len = 1, pattern = "^[^._]+$")
  assertCharacter(x = layer, null.ok = TRUE)

  if(is.null(layer)) layer <- names(obj@layers)[1]
  if(!(layer %in% names(obj@layers))){
    stop("layer '", layer, "' not found in 'obj'.")
  }
  if(!is.null(cost) && !(cost %in% names(obj@layers))){
    stop("cost layer '", cost, "' not found in 'obj'.")
  }
  # msr() reads a name that matches a layer as that layer
  if(name %in% names(obj@layers)){
    stop("'name' is '", name, "', which is a layer of 'obj'; choose another.",
         call. = FALSE)
  }

  theRes <- msk_res(obj)
  # a distance in cells becomes one in map units only if a cell is as wide as
  # it is high
  if(unit == "map" && theRes[1] != theRes[2]){
    stop("unit = \"map\" needs square cells, but the resolution is ",
         theRes[1], " x ", theRes[2], ".", call. = FALSE)
  }

  dims <- obj@dims
  vals <- msk_pull(obj, layer)
  surface <- if(is.null(cost)) NULL else msk_pull(obj, cost)

  uVals <- sort(unique(vals[!is.na(vals)]))
  n <- length(uVals)
  cells <- lapply(uVals, function(k) which(!is.na(vals) & vals == k))

  # the diagonal stays Inf: a class has no distance to itself, and this keeps
  # min(distance.others) reading as the nearest neighbour; unreachable pairs
  # become NA
  mat <- matrix(Inf, nrow = n, ncol = n)

  for(p in seq_len(n)){

    if(routing == "straight"){
      # the transform gives every cell its distance to the nearest cell of
      # class p, so the nearest cell of class q is where it is smallest
      src <- rep(0, length(vals))
      src[cells[[p]]] <- 1
      d <- sqrt(distanceCpp(vals = src, nrow = dims[2], ncol = dims[1],
                            method = "euclidean"))
      for(q in seq_len(n)){
        if(p == q) next
        if(is.null(cost)){
          mat[p, q] <- min(d[cells[[q]]])
        } else {
          # the closest pair: that cell of q and the cell of p nearest to it
          to <- cells[[q]][which.min(d[cells[[q]]])]
          mat[p, q] <- .straight_path(cells[[p]], to, surface, dims)
        }
      }
      next
    }

    cellCost <- if(is.null(cost)) rep(1, length(vals)) else surface

    # the accumulated cost of the cheapest path, entering the target class at
    # its cheapest cell
    cd <- costDistanceCpp(cost = cellCost, from = cells[[p]],
                          nrow = dims[2], ncol = dims[1], diagonal = TRUE)
    for(q in seq_len(n)){
      if(p == q) next
      reached <- cd$dist[cells[[q]]]
      reached <- reached[is.finite(reached)]
      mat[p, q] <- if(length(reached) == 0) NA_real_ else min(reached)
    }
  }

  if(unit == "map") mat <- mat * theRes[1]

  obj <- .store_class(obj, layer, name, as.integer(uVals), mat)

  obj <- .update_mosaik(obj, step = step)

  return(obj)
}

#' Sum the costs along the direct line between two classes
#'
#' Walks the straight line from the cell of the source class nearest to the
#' target cell, which is the target class's cell nearest to the source class,
#' and sums the costs of the cells it crosses, whether or not a cheaper detour
#' exists. Only one target cell is compared with the source cells, so memory
#' grows with the size of the classes, not with their product.
#'
#' @param fromCells [`integer(.)`][integer]\cr cell indices of the source class.
#' @param toCell [`integer(1)`][integer]\cr the target class's cell nearest to
#'   the source class.
#' @param cellCost [`numeric(.)`][numeric]\cr the cost of each cell.
#' @param dims [`integer(2)`][integer]\cr grid dimensions, columns then rows.
#' @return A single numeric value, \code{NA} when the line crosses a cell that
#'   cannot be crossed.
#' @keywords internal

.straight_path <- function(fromCells, toCell, cellCost, dims){

  ncol <- dims[1]
  fr <- (fromCells - 1) %/% ncol; fc <- (fromCells - 1) %% ncol
  r1 <- (toCell - 1) %/% ncol;    c1 <- (toCell - 1) %% ncol

  # the source cell nearest to the target cell: with the target cell nearest
  # to the source class, the two are a closest pair
  best <- which.min((fr - r1)^2 + (fc - c1)^2)
  r0 <- fr[best]; c0 <- fc[best]

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
