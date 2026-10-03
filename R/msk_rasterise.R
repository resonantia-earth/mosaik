#' Rasterise points, lines or polygons
#'
#' Write vector geometries into a layer: each cell a geometry covers gets the
#' geometry's \code{id}, all other cells \code{NA}.
#'
#' @param obj [`mosaik`]\cr the mosaik whose grid the geometries are written
#'   on; an empty grid from \code{\link{mosaik}(extent, res)} will do.
#' @param geom [`data.frame`][data.frame]\cr the geometries, one vertex per row,
#'   with the numeric columns \code{x} and \code{y} (the coordinates, in the
#'   units of the mosaik's extent), \code{id} (the geometry the vertex belongs
#'   to) and, optionally, \code{part} (the part of a multi-geometry the vertex
#'   belongs to). \code{\link{msk_spaghettify}} makes this table from sf or
#'   terra vectors.
#' @param type [`character(1)`][character]\cr what the geometries are:
#'   \code{"polygon"} (default), \code{"line"} or \code{"point"}.
#' @param name [`character(1)`][character]\cr the layer name. Default
#'   \code{"values"}.
#' @return A mosaik in which each cell carries the id of the geometry covering
#'   it.
#' @details
#'   The rows of one \code{id} form one geometry, and all rows of one
#'   \code{part} within it form one piece of that geometry, in the order they
#'   are given. Without a \code{part} column, every \code{id} is one piece. A
#'   multipolygon or multiline is therefore one \code{id} with several parts;
#'   to treat its pieces as separate geometries, give each its own \code{id}
#'   beforehand.
#'
#'   \describe{
#'     \item{\code{"point"}}{a point covers the cell it lies in.}
#'     \item{\code{"line"}}{a line covers the cells it passes through,
#'       between consecutive vertices of a part.}
#'     \item{\code{"polygon"}}{each part is a ring, closed automatically if its
#'       last vertex is not its first. A polygon covers the cells whose centre
#'       lies inside an odd number of its rings, so a ring inside another one is
#'       a hole, and a ring inside a hole is an island again. Holes need no
#'       marking of their own.}
#'   }
#'
#'   Where geometries overlap, the one given later wins. Values of the
#'   geometries, such as an attribute of each polygon, are written onto the
#'   cells afterwards with \code{\link{mdf_replace}}\code{(old = id, new =
#'   value)}.
#' @family utilities
#' @examples
#' m <- mosaik(extent = c(0, 60, 0, 60), res = 1)
#'
#' # two polygons, the second with a hole and an island in it
#' polygons <- data.frame(
#'   x    = c(5, 25, 25, 5,   30, 55, 55, 30,   36, 49, 49, 36,   40, 45, 45, 40),
#'   y    = c(5, 5, 25, 25,   30, 30, 55, 55,   36, 36, 49, 49,   40, 40, 45, 45),
#'   id   = rep(c(1, 2), c(4, 12)),
#'   part = rep(c(1, 1, 2, 3), each = 4))
#'
#' # a line in two parts, a straight one and a curve, and three points
#' t <- seq(0, 1, length.out = 25)
#' curve_x <- 10 + 45 * t
#' curve_y <- 5 + 20 * sin(pi * t) * (1 - 0.4 * t)
#' lines <- data.frame(x = c(5, 50, curve_x), y = c(55, 40, curve_y), id = 1,
#'                     part = rep(c(1, 2), c(2, length(t))))
#' points <- data.frame(x = c(10, 30, 50), y = c(45, 15, 5), id = 1:3)
#'
#' m <- m |>
#'   msk_rasterise(geom = polygons, name = "polygons") |>
#'   msk_rasterise(geom = lines, type = "line", name = "lines") |>
#'   msk_rasterise(geom = points, type = "point", name = "points")
#' msk_vis(m, .layer("polygons"), .layer("lines"), .layer("points"))
#' @importFrom checkmate assertClass assertDataFrame assertNames assertNumeric
#'   assertChoice assertCharacter
#' @export

msk_rasterise <- function(obj,
                          geom,
                          type = "polygon",
                          name = "values") {

  step <- .step()

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertDataFrame(x = geom, min.rows = 1)
  assertNames(x = names(geom), must.include = c("x", "y", "id"))
  assertNumeric(x = geom$x, any.missing = FALSE)
  assertNumeric(x = geom$y, any.missing = FALSE)
  assertNumeric(x = geom$id, any.missing = FALSE)
  if (!is.null(geom$part)) assertNumeric(x = geom$part, any.missing = FALSE)
  assertChoice(x = type, choices = c("point", "line", "polygon"))
  assertCharacter(x = name, len = 1, any.missing = FALSE)

  # pull data ----
  dims <- obj@dims
  ext  <- obj@extent
  res  <- msk_res(obj)
  if (is.null(geom$part)) geom$part <- 1

  # the pieces of each geometry, in the order they are given
  ids <- unique(geom$id)
  piece <- paste(geom$id, geom$part)
  pieces <- split(seq_len(nrow(geom)), factor(piece, levels = unique(piece)))
  piece_id <- vapply(pieces, function(r) as.numeric(geom$id[r[1]]), numeric(1))

  temp <- rep(NA_real_, prod(dims))

  # body ----
  if (type == "point") {
    cells <- .cell(obj, x = geom$x, y = geom$y)
    keep <- !is.na(cells)
    temp[cells[keep]] <- geom$id[keep]

  } else if (type == "line") {
    # sample each segment at a tenth of a cell, so it meets every cell it
    # passes through except where it only grazes a corner
    step_len <- min(res) / 10
    for (rows in pieces) {
      if (length(rows) < 2) {
        stop("a line needs at least two vertices; id ", geom$id[rows[1]],
             " has a part with one.", call. = FALSE)
      }
      x <- geom$x[rows]; y <- geom$y[rows]
      for (k in seq_len(length(rows) - 1)) {
        n <- max(2, ceiling(sqrt((x[k + 1] - x[k])^2 + (y[k + 1] - y[k])^2) /
                              step_len) + 1)
        t <- seq(0, 1, length.out = n)
        cells <- .cell(obj, x = x[k] + t * (x[k + 1] - x[k]),
                          y = y[k] + t * (y[k + 1] - y[k]))
        temp[cells[!is.na(cells)]] <- geom$id[rows[1]]
      }
    }

  } else {
    # cell centres, row by row from the top-left
    xs <- seq(ext[1] + res[1] / 2, by = res[1], length.out = dims[1])
    ys <- seq(ext[4] - res[2] / 2, by = -res[2], length.out = dims[2])
    pts <- cbind(rep(xs, times = dims[2]), rep(ys, each = dims[1]))

    for (i in ids) {
      # even-odd across the rings of one polygon: a hole is a ring inside
      # another, so no ring needs to be marked as one
      crossed <- integer(prod(dims))
      for (rows in pieces[piece_id == i]) {
        ring <- cbind(geom$x[rows], geom$y[rows])
        if (any(ring[1, ] != ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
        if (nrow(ring) < 4) {
          stop("a polygon ring needs at least three vertices; id ", i,
               " has a part with fewer.", call. = FALSE)
        }
        inside <- pointInPolyCpp(vert = pts, geom = ring, invert = FALSE)
        crossed <- crossed + (inside >= 1L)
      }
      temp[crossed %% 2 == 1] <- i
    }
  }

  # build output ----
  .update_mosaik(obj, values = temp, keep = FALSE, step = step)
}
