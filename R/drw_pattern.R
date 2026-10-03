#' Make a geometric pattern surface
#'
#' Generates regular geometric patterns such as checkerboards, stripes,
#' concentric rings, sinusoidal waves, block grids, or hexagonal tilings.
#'
#' @param obj [`mosaik`]\cr the mosaik whose grid the pattern is generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [`character(1)`][character]\cr the pattern type: \code{"checkerboard"}
#'   (default), \code{"stripes"}, \code{"rings"}, \code{"waves"}, \code{"grid"},
#'   or \code{"hexagonal"}.
#' @param frequency [`numeric(1)`][numeric]\cr how fine the pattern is. For the
#'   tiles (\code{"checkerboard"}, \code{"grid"}, \code{"hexagonal"}), the tile
#'   size in cells; for hexagons their width across the flat sides, which gives
#'   exactly symmetric hexagons when it is even. For \code{"stripes"},
#'   \code{"rings"} and \code{"waves"},
#'   the number of stripes, rings or waves across the grid (from the centre to
#'   the corners, for rings).
#' @param n [`integerish(1)`][integer]\cr number of distinct classes for tiled
#'   patterns (\code{"checkerboard"}, \code{"hexagonal"}).
#'   Defaults to the minimum proper colouring: 2 for checkerboard, 3 for
#'   hexagonal. Set to \code{0} to give each tile a unique sequential ID
#'   (useful for zoning). Ignored by other pattern types.
#' @param angle [`numeric(1)`][numeric]\cr the direction in which stripes and
#'   waves repeat, in degrees counter-clockwise from the x-axis (default 0:
#'   left to right, so the bands run vertically). Applies to \code{"stripes"}
#'   and \code{"waves"}.
#' @param name [`character(1)`][character]\cr the layer name. Default
#'   \code{"values"}.
#' @return A \code{\link{mosaik}} with pattern values.
#' @family generator functions
#' @examples
#' m <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # tiles of 10 cells: checkerboard, grid and hexagons
#' p <- m |>
#'   drw_pattern(type = "checkerboard", frequency = 10, name = "checkerboard") |>
#'   drw_pattern(type = "grid", frequency = 10, name = "grid") |>
#'   drw_pattern(type = "hexagonal", frequency = 10, name = "hexagonal")
#' msk_vis(p, .layer("checkerboard"), .layer("grid"), .layer("hexagonal"))
#'
#' # each tile numbered on its own (n = 0), e.g. as zones; drawn in shuffled
#' # colours, since neighbouring numbers would otherwise get similar ones
#' p <- m |>
#'   drw_pattern(type = "checkerboard", frequency = 20, n = 0,
#'               name = "checkerboard") |>
#'   drw_pattern(type = "hexagonal", frequency = 10, n = 0, name = "hexagonal")
#' set.seed(1)
#' shuffled <- function(layer) {
#'   sample(hcl.colors(max(msk_pull(p, layer)), "Dark 3"))
#' }
#' msk_vis(p,
#'         .layer("checkerboard", colours = shuffled("checkerboard"),
#'                legend = FALSE),
#'         .layer("hexagonal", colours = shuffled("hexagonal"), legend = FALSE))
#'
#' # 5 stripes, rings and waves across the grid, stripes and waves at 60 degrees
#' p <- m |>
#'   drw_pattern(type = "stripes", frequency = 5, angle = 60, name = "stripes") |>
#'   drw_pattern(type = "rings", frequency = 5, name = "rings") |>
#'   drw_pattern(type = "waves", frequency = 5, angle = 60, name = "waves")
#' msk_vis(p, .layer("stripes"), .layer("rings"), .layer("waves"))
#' @importFrom checkmate assertClass assertIntegerish assertChoice assertNumber
#'   assertCharacter
#' @export

drw_pattern <- function(obj, type = "checkerboard", frequency = 10,
                         n = NULL, angle = 0, name = "values") {

  step <- .step()

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type,
               choices = c("checkerboard", "stripes", "rings", "waves",
                           "grid", "hexagonal"))
  assertNumber(x = frequency, lower = 1)
  assertIntegerish(x = n, len = 1, lower = 0, null.ok = TRUE)
  assertNumber(x = angle, finite = TRUE)
  assertCharacter(x = name, len = 1)

  if (is.null(n)) {
    n <- switch(type, checkerboard = 2L, hexagonal = 3L, 0L)
  } else {
    n <- as.integer(n)
  }

  ncols <- obj@dims[1]
  nrows <- obj@dims[2]

  col_idx <- rep(seq_len(ncols), times = nrows) - 1L
  row_idx <- rep(rev(seq_len(nrows)), each = ncols) - 1L

  vals <- switch(type,

    checkerboard = {
      block <- as.integer(frequency)
      tile_q <- floor(col_idx / block)
      tile_r <- floor(row_idx / block)
      if (n == 0L) {
        n_tiles_x <- ceiling(ncols / block)
        as.numeric(tile_r * n_tiles_x + tile_q + 1L)
      } else {
        as.numeric((tile_q + tile_r) %% n)
      }
    },

    stripes = {
      # 'frequency' stripes of 1 across the grid, each followed by a band of 0
      # of the same width
      rad <- angle * pi / 180
      proj <- col_idx * cos(rad) + row_idx * sin(rad)
      phase <- (proj - min(proj)) / (diff(range(proj)) + 1) * frequency
      as.numeric(floor(phase * 2) %% 2 == 0)
    },

    rings = {
      cx <- (ncols - 1) / 2
      cy <- (nrows - 1) / 2
      dist <- sqrt((col_idx - cx)^2 + (row_idx - cy)^2)
      max_dist <- sqrt(cx^2 + cy^2)
      scaled <- dist / max_dist * frequency
      (cos(scaled * 2 * pi) + 1) / 2
    },

    waves = {
      rad <- angle * pi / 180
      proj <- col_idx * cos(rad) + row_idx * sin(rad)
      max_proj <- max(abs(proj))
      if (max_proj == 0) max_proj <- 1
      scaled <- proj / max_proj * frequency
      (sin(scaled * 2 * pi) + 1) / 2
    },

    grid = {
      block <- as.integer(frequency)
      x_edge <- (col_idx %% block) == 0 | col_idx == (ncols - 1L)
      y_edge <- (row_idx %% block) == 0 | row_idx == (nrows - 1L)
      as.numeric(x_edge | y_edge)
    },

    hexagonal = {
      # pointy-topped hexagons 'frequency' cells wide across their flat sides.
      # Measured from the cell centres, the hexagon centres fall on cell
      # corners, so with an even width every boundary between two hexagons of
      # a row runs along cell edges and each hexagon is symmetric about its
      # vertical axis
      w <- as.numeric(frequency)
      s <- w / sqrt(3)
      x <- col_idx + 0.5
      y <- row_idx + 0.5
      q <- (x * sqrt(3) / 3 - y / 3) / s
      r_hex <- (2 / 3 * y) / s

      rx <- round(q)
      rz <- round(r_hex)
      ry <- round(-q - r_hex)

      xd <- abs(rx - q)
      yd <- abs(ry - (-q - r_hex))
      zd <- abs(rz - r_hex)

      fix_x <- (xd > yd) & (xd > zd)
      fix_z <- (!fix_x) & (zd > yd)

      rx <- ifelse(fix_x, -ry - rz, rx)
      rz <- ifelse(fix_z, -rx - ry, rz)

      if (n == 0L) {
        ids <- paste0(rx, ",", rz)
        as.numeric(match(ids, unique(ids)))
      } else {
        as.numeric(((2L * rx + rz) %% n + n) %% n)
      }
    }
  )


  .update_mosaik(obj, values = vals, step = step)
}
