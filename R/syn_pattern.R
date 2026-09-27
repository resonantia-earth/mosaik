# TODO: the hexagonal pattern doesn't look particularly good. It's not regular,
# per se, the hexagons are not the same dimensions in all directions, for
# instance, they are 8 cells to the right and 9 to the left from the middle axis.

#' Make a geometric pattern surface
#'
#' Generates regular geometric patterns such as checkerboards, stripes,
#' concentric rings, sinusoidal waves, block grids, or hexagonal tilings.
#'
#' @param obj [mosaik]\cr the mosaik whose grid the pattern is generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [character(1)][character]\cr the pattern type: \code{"checkerboard"}
#'   (default), \code{"stripes"}, \code{"rings"}, \code{"waves"}, \code{"grid"},
#'   or \code{"hexagonal"}.
#' @param frequency [numeric(1)][numeric]\cr controls the spatial frequency of
#'   the pattern. For \code{"checkerboard"} and \code{"grid"}, the block size
#'   in cells. For \code{"stripes"}, the stripe width. For \code{"rings"}, the
#'   number of rings. For \code{"waves"}, the number of wave cycles across the
#'   grid. For \code{"hexagonal"}, the tile size in cells.
#' @param n [integerish(1)][integer]\cr number of distinct classes for tiled
#'   patterns (\code{"checkerboard"}, \code{"hexagonal"}).
#'   Defaults to the minimum proper colouring: 2 for checkerboard, 3 for
#'   hexagonal. Set to \code{0} to give each tile a unique sequential ID
#'   (useful for zoning). Ignored by other pattern types.
#' @param angle [numeric(1)][numeric]\cr rotation angle in degrees (default 0).
#'   Applies to \code{"stripes"} and \code{"waves"}.
#' @param name [character(1)][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param seed [integerish(1)][integer]\cr ignored (patterns are deterministic),
#'   accepted for API consistency.
#' @return A \code{\link{mosaik}} with pattern values.
#' @family generator functions
#' @examples
#' m <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # checkerboard with 10-cell blocks
#' syn_pattern(m, type = "checkerboard", frequency = 10)
#'
#' # diagonal stripes
#' syn_pattern(m, type = "stripes", frequency = 8, angle = 60)
#'
#' # hexagonal tiling
#' syn_pattern(m, type = "hexagonal", frequency = 10)
#' @importFrom checkmate assertClass assertIntegerish assertChoice assertNumber
#'   assertCharacter
#' @export

syn_pattern <- function(obj, type = "checkerboard", frequency = 10,
                         n = NULL, angle = 0, name = "values",
                         seed = NULL) {

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type,
               choices = c("checkerboard", "stripes", "rings", "waves",
                           "grid", "hexagonal"))
  assertNumber(x = frequency, lower = 1)
  assertIntegerish(x = n, len = 1, lower = 0, null.ok = TRUE)
  assertNumber(x = angle, finite = TRUE)
  assertCharacter(x = name, len = 1)
  assertIntegerish(x = seed, len = 1, null.ok = TRUE)

  if (is.null(n)) {
    n <- switch(type, checkerboard = 2L, hexagonal = 3L, 0L)
  } else {
    n <- as.integer(n)
  }

  ncols <- obj@dims[1]
  nrows <- obj@dims[2]

  col_idx <- rep(seq_len(ncols), times = nrows) - 1L
  row_idx <- rep(seq_len(nrows), each = ncols) - 1L

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
      rad <- angle * pi / 180
      proj <- col_idx * cos(rad) + row_idx * sin(rad)
      width <- as.numeric(frequency)
      as.numeric(floor(proj / width) %% 2)
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
      s <- as.numeric(frequency)
      q <- (col_idx * sqrt(3) / 3 - row_idx / 3) / s
      r_hex <- row_idx / (1.5 * s)

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

  prov <- msk_prov("syn_pattern",
                     list(type = type, frequency = frequency,
                          n = n, angle = angle, name = name))

  msk_set(obj, name, vals, prov)
}
