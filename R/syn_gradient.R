# TODO: not sure it makes sense to have both obj and origin. I think they could
# be merged because what would happen if obj is different than oigin? Why would
# we allow this?
# Why does it have seed when that is anyway ignored?
# I see that position starts at the top left, while the coordinate system at the
# lower left. I think both should be at the lower left.
# could "angle" not also make sense for the shapes?

#' Make a gradient surface
#'
#' A gradient is the distance from an origin, scaled to \[0, 1\]. The origin
#' can be a point, a line, a geometric shape, or an arbitrary binary mosaik.
#'
#' @param obj [mosaik]\cr the mosaik whose grid the gradient is generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param origin [mosaik]\cr an optional binary mosaik whose foreground cells
#'   serve as the origin. Overrides \code{type} if given.
#' @param type [character(1)][character]\cr the type of gradient:
#'   \code{"planar"} (default), \code{"point"}, \code{"line"},
#'   \code{"circle"}, \code{"rectangle"}, \code{"square"}, \code{"polygon"},
#'   \code{"ellipse"}, \code{"triangle"}, \code{"hexagon"}.
#'   Ignored when \code{origin} is provided.
#' @param angle [numeric(1)][numeric]\cr rotation angle in degrees. For
#'   \code{"planar"}, controls the gradient direction (0 = left-to-right,
#'   90 = bottom-to-top). For \code{"line"}, the line angle.
#' @param position [numeric(2)][numeric]\cr the centre position of the origin
#'   shape as \code{c(x, y)} in grid-relative coordinates (0-1). Default
#'   \code{c(0.5, 0.5)} (centre). Used by \code{"point"}, \code{"circle"},
#'   \code{"rectangle"}, \code{"square"}, \code{"ellipse"}, \code{"triangle"},
#'   \code{"hexagon"}.
#' @param size [numeric(1)][numeric]\cr size of the origin shape as a fraction
#'   of the shorter grid dimension (0-1). Default 0.3. Used by all shape types
#'   except \code{"planar"}, \code{"point"}, and \code{"line"}.
#' @param invert [logical(1)][logical]\cr if \code{TRUE}, invert the gradient
#'   so that cells near the origin have high values and distant cells have low
#'   values. Default \code{FALSE}.
#' @param name [character(1)][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param seed [integerish(1)][integer]\cr ignored (gradients are deterministic).
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
#' # left-to-right planar gradient (default)
#' syn_gradient(m)
#'
#' # diagonal gradient
#' syn_gradient(m, type = "planar", angle = 45)
#'
#' # radial gradient from a point
#' syn_gradient(m, type = "point", position = c(0.3, 0.7))
#'
#' # distance from a shape
#' syn_gradient(m, type = "circle", size = 0.2)
#'
#' # inverted: high values near the origin
#' syn_gradient(m, type = "point", invert = TRUE)
#'
#' @importFrom checkmate assertIntegerish assertChoice assertNumber
#'   assertLogical assertCharacter assertClass
#' @export

syn_gradient <- function(obj, origin = NULL, type = "planar",
                          angle = 0, position = c(0.5, 0.5), size = 0.3,
                          invert = FALSE, name = "values",
                          seed = NULL) {

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type,
               choices = c("planar", "point", "line", "rectangle", "square",
                           "polygon", "ellipse", "circle",
                           "triangle", "hexagon"))
  assertNumber(x = angle, finite = TRUE)
  assertNumeric(x = position, len = 2, lower = 0, upper = 1)
  assertNumber(x = size, lower = 0, upper = 1)
  assertLogical(x = invert, len = 1)
  assertCharacter(x = name, len = 1)
  assertIntegerish(x = seed, len = 1, null.ok = TRUE)

  ncols <- obj@dims[1]
  nrows <- obj@dims[2]
  n_cells <- ncols * nrows

  if (!is.null(origin)) {

    assertClass(x = origin, classes = "mosaik")
    binary <- msk_pull(origin)
    if (!isBinaryCpp(vals = binary)) {
      stop("'origin' must be a binary mosaik (values 0 and 1).")
    }
    temp <- distanceCpp(vals = binary, nrow = nrows, ncol = ncols,
                        method = "euclidean")
    temp <- sqrt(temp)

  } else if (type == "planar") {

    rad <- angle * pi / 180
    col_idx <- rep(seq_len(ncols), times = nrows) - 1L
    row_idx <- rep(seq_len(nrows), each = ncols) - 1L
    temp <- col_idx * cos(rad) + row_idx * sin(rad)

  } else if (type == "point") {

    cx <- position[1] * ncols
    cy <- position[2] * nrows
    col_idx <- rep(seq_len(ncols), times = nrows) - 0.5
    row_idx <- rep(seq_len(nrows), each = ncols) - 0.5
    temp <- sqrt((col_idx - cx)^2 + (row_idx - cy)^2)

  } else if (type == "line") {

    rad <- angle * pi / 180
    cx <- position[1] * ncols
    cy <- position[2] * nrows
    nx <- -sin(rad)
    ny <- cos(rad)
    col_idx <- rep(seq_len(ncols), times = nrows) - 0.5
    row_idx <- rep(seq_len(nrows), each = ncols) - 0.5
    temp <- abs((col_idx - cx) * nx + (row_idx - cy) * ny)

  } else {

    # shape-based origins: rasterise the shape as a binary mask, then distance
    cx <- position[1] * ncols
    cy <- position[2] * nrows
    short <- min(ncols, nrows)
    radius <- size * short / 2

    col_idx <- rep(seq_len(ncols), times = nrows) - 0.5
    row_idx <- rep(seq_len(nrows), each = ncols) - 0.5

    binary <- switch(type,
      circle = {
        dist <- sqrt((col_idx - cx)^2 + (row_idx - cy)^2)
        as.numeric(dist <= radius)
      },
      rectangle = {
        hw <- radius * 1.5
        hh <- radius
        as.numeric(abs(col_idx - cx) <= hw & abs(row_idx - cy) <= hh)
      },
      square = {
        as.numeric(abs(col_idx - cx) <= radius & abs(row_idx - cy) <= radius)
      },
      ellipse = {
        a <- radius * 1.5
        b <- radius
        dist <- ((col_idx - cx) / a)^2 + ((row_idx - cy) / b)^2
        as.numeric(dist <= 1)
      },
      triangle = {
        h <- radius * sqrt(3)
        ry <- row_idx - cy
        rx <- col_idx - cx
        inside <- (ry >= -h / 3) &
                  (ry <= 2 * h / 3 - abs(rx) * sqrt(3))
        as.numeric(inside)
      },
      hexagon = {
        rx <- abs(col_idx - cx)
        ry <- abs(row_idx - cy)
        inside <- (rx <= radius) &
                  (ry <= radius * sqrt(3) / 2) &
                  (radius * sqrt(3) - sqrt(3) * rx - ry >= 0)
        as.numeric(inside)
      },
      polygon = {
        rx <- abs(col_idx - cx)
        ry <- abs(row_idx - cy)
        cut <- radius * cos(pi / 8)
        inside <- (rx <= cut) & (ry <= cut) &
                  (rx + ry <= radius * (1 + cos(pi / 8)))
        as.numeric(inside)
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

  prov <- msk_prov("syn_gradient",
                     list(type = type, angle = angle, position = position,
                          size = size, invert = invert, name = name))

  msk_set(obj, name, temp, prov)
}
