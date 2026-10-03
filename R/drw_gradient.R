#' Make a gradient surface
#'
#' A gradient is the distance from an origin, scaled to \[0, 1\]. The origin
#' can be a plane, a point, a line, a circle, a square, or any binary layer of
#' the mosaik, which covers every other shape.
#'
#' @param obj [`mosaik`]\cr the mosaik whose grid the gradient is generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [`character(1)`][character]\cr the type of gradient:
#'   \code{"planar"} (default), \code{"point"}, \code{"line"},
#'   \code{"circle"} or \code{"square"}. Ignored when \code{origin} is given.
#' @param angle [`numeric(1)`][numeric]\cr angle in degrees, counter-clockwise
#'   from the x-axis. For \code{"planar"}, the direction in which the gradient
#'   increases (0 = left to right, 90 = bottom to top). For \code{"line"}, the
#'   direction of the line. For \code{"square"}, the turn of the square around
#'   its centre.
#' @param position [`numeric(2)`][numeric]\cr the centre of the origin as
#'   \code{c(x, y)} in grid-relative coordinates (0-1), counted from the
#'   lower-left corner like the map coordinates. Default
#'   \code{c(0.5, 0.5)} (centre). Used by \code{"point"}, \code{"line"},
#'   \code{"circle"} and \code{"square"}.
#' @param size [`numeric(1)`][numeric]\cr size of the circle (its diameter) or
#'   the square (its side) as a fraction of the shorter grid dimension (0-1).
#'   Default 0.3.
#' @param origin [`character(1)`][character]\cr the name of a binary layer of
#'   \code{obj} whose cells of value 1 are the origin, for any shape other than
#'   those of \code{type}. Overrides \code{type} if given.
#' @param invert [`logical(1)`][logical]\cr if \code{TRUE}, invert the gradient
#'   so that cells near the origin have high values and distant cells have low
#'   values. Default \code{FALSE}.
#' @param name [`character(1)`][character]\cr the layer name. Default
#'   \code{"values"}.
#' @return A \code{\link{mosaik}} with gradient values scaled between 0 and 1.
#' @references
#'   Meijster, A., Roerdink, J.B.T.M., Hesselink, W.H., 2000. A general
#'   algorithm for computing distance transforms in linear time, in: Goutsias,
#'   J., Vincent, L., Bloomberg, D.S. (Eds.), Mathematical Morphology and Its
#'   Applications to Image and Signal Processing. Springer, pp. 331-340.
#' @family generator functions
#' @examples
#' m <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # a plane, left to right (default) and diagonal
#' g <- m |>
#'   drw_gradient(name = "planar") |>
#'   drw_gradient(angle = 45, name = "planar_45")
#' msk_vis(g, .layer("planar"), .layer("planar_45"))
#'
#' # from a point, and from a line at 30 degrees through the centre
#' g <- m |>
#'   drw_gradient(type = "point", position = c(0.3, 0.7), name = "point") |>
#'   drw_gradient(type = "line", angle = 30, name = "line")
#' msk_vis(g, .layer("point"), .layer("line"))
#'
#' # from a circle, a square, and a square turned by 30 degrees
#' g <- m |>
#'   drw_gradient(type = "circle", size = 0.4, name = "circle") |>
#'   drw_gradient(type = "square", size = 0.4, name = "square") |>
#'   drw_gradient(type = "square", size = 0.4, angle = 30,
#'                name = "square_30")
#' msk_vis(g, .layer("circle"), .layer("square"), .layer("square_30"))
#'
#' # inverted: high values near the origin
#' g <- drw_gradient(m, type = "point", invert = TRUE, name = "inverted")
#' msk_vis(g, .layer("inverted"))
#'
#' # any other shape as the origin: here the forest of the example data
#' g <- landscape |>
#'   mdf_binarise(match = 47, layer = "cover", add = "forest") |>
#'   drw_gradient(origin = "forest", name = "from_forest")
#' msk_vis(g, .layer("forest"), .layer("from_forest"))
#' @importFrom checkmate assertIntegerish assertChoice assertNumber
#'   assertLogical assertCharacter assertClass
#' @export

drw_gradient <- function(obj, type = "planar", angle = 0,
                          position = c(0.5, 0.5), size = 0.3, origin = NULL,
                          invert = FALSE, name = "values") {

  step <- .step()

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type,
               choices = c("planar", "point", "line", "circle", "square"))
  assertNumber(x = angle, finite = TRUE)
  assertNumeric(x = position, len = 2, lower = 0, upper = 1)
  assertNumber(x = size, lower = 0, upper = 1)
  assertCharacter(x = origin, len = 1, null.ok = TRUE)
  assertLogical(x = invert, len = 1)
  assertCharacter(x = name, len = 1)

  ncols <- obj@dims[1]
  nrows <- obj@dims[2]
  n_cells <- ncols * nrows

  if (!is.null(origin)) {

    if (!origin %in% names(obj@layers)) {
      stop("'origin' names no layer of 'obj': ", origin, call. = FALSE)
    }
    binary <- msk_pull(obj, origin)
    if (!isBinaryCpp(vals = binary)) {
      stop("'origin' must be a binary layer (values 0 and 1).", call. = FALSE)
    }
    temp <- distanceCpp(vals = binary, nrow = nrows, ncol = ncols,
                        method = "euclidean")
    temp <- sqrt(temp)

  } else if (type == "planar") {

    rad <- angle * pi / 180
    col_idx <- rep(seq_len(ncols), times = nrows) - 1L
    row_idx <- rep(rev(seq_len(nrows)), each = ncols) - 1L
    temp <- col_idx * cos(rad) + row_idx * sin(rad)

  } else if (type == "point") {

    cx <- position[1] * ncols
    cy <- position[2] * nrows
    col_idx <- rep(seq_len(ncols), times = nrows) - 0.5
    row_idx <- rep(rev(seq_len(nrows)), each = ncols) - 0.5
    temp <- sqrt((col_idx - cx)^2 + (row_idx - cy)^2)

  } else if (type == "line") {

    rad <- angle * pi / 180
    cx <- position[1] * ncols
    cy <- position[2] * nrows
    nx <- -sin(rad)
    ny <- cos(rad)
    col_idx <- rep(seq_len(ncols), times = nrows) - 0.5
    row_idx <- rep(rev(seq_len(nrows)), each = ncols) - 0.5
    temp <- abs((col_idx - cx) * nx + (row_idx - cy) * ny)

  } else {

    # shape-based origins: rasterise the shape as a binary mask, then distance
    cx <- position[1] * ncols
    cy <- position[2] * nrows
    short <- min(ncols, nrows)
    radius <- size * short / 2

    col_idx <- rep(seq_len(ncols), times = nrows) - 0.5
    row_idx <- rep(rev(seq_len(nrows)), each = ncols) - 0.5

    binary <- switch(type,
      circle = {
        dist <- sqrt((col_idx - cx)^2 + (row_idx - cy)^2)
        as.numeric(dist <= radius)
      },
      square = {
        # turned by 'angle' around its centre
        rad <- angle * pi / 180
        u <- (col_idx - cx) * cos(rad) + (row_idx - cy) * sin(rad)
        v <- -(col_idx - cx) * sin(rad) + (row_idx - cy) * cos(rad)
        as.numeric(abs(u) <= radius & abs(v) <= radius)
      }
    )

    temp <- distanceCpp(vals = binary, nrow = nrows, ncol = ncols,
                        method = "euclidean")
    temp <- sqrt(temp)
  }

  # scale to [0, 1]
  rng <- range(temp, na.rm = TRUE)
  if (rng[2] > rng[1]) {
    temp <- (temp - rng[1]) / (rng[2] - rng[1])
  } else {
    temp <- rep(0, n_cells)
  }

  if (invert) temp <- 1 - temp


  .update_mosaik(obj, values = temp, step = step)
}
