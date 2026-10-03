#' Make a tessellation surface
#'
#' Partition the grid into discrete regions. \code{"voronoi"} places random seed
#' points and assigns each cell to the nearest seed. \code{"gibbs"} spaces the
#' seeds more regularly before the same assignment: a seed closer than the
#' inhibition radius to an earlier one is rejected (a hard-core process, the
#' Strauss process with complete inhibition). \code{"rectangle"} is the random rectangular cluster: rectangles
#' of random size are dropped at random places, each covering the ones below,
#' until the grid is covered (a "dead leaves" model).
#'
#' @param obj [`mosaik`]\cr the mosaik whose grid the tessellation is generated
#'   on. Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [`character(1)`][character]\cr tessellation type: \code{"voronoi"}
#'   (default), \code{"gibbs"} or \code{"rectangle"}.
#' @param n [`integerish(1)`][integer]\cr for \code{"voronoi"} and
#'   \code{"gibbs"}: the number of seed points, and so of regions. Default 20.
#' @param interaction [`numeric(1)`][numeric]\cr for \code{"gibbs"} only: the
#'   inhibition radius as a fraction of the grid diagonal; no two seeds are
#'   closer than this. If not all \code{n} seeds fit, a warning says how many
#'   were placed. Default 0.1.
#' @param size [`integerish(2)`][integer]\cr for \code{"rectangle"} only: the
#'   smallest and largest side length of a rectangle, in cells. Each side is
#'   drawn at random from this range. Default \code{c(5, 20)}.
#' @param name [`character(1)`][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param seed [`integerish(1)`][integer]\cr random seed for reproducibility.
#' @return A \code{\link{mosaik}} with integer region labels.
#' @references
#' Aurenhammer F. Voronoi diagrams — a survey of a fundamental geometric data
#' structure. ACM Computing Surveys. 1991;23(3):345-405.
#' \href{https://doi.org/10.1145/116873.116880}{10.1145/116873.116880}
#'
#' Strauss DJ. A model for clustering. Biometrika. 1975;62(2):467-475.
#'
#' Ripley BD. Modelling spatial patterns (with discussion). Journal of the Royal
#' Statistical Society, Series B. 1977;39(2):172-212.
#' @family generator functions
#' @examples
#' m <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # Voronoi cells around 20 random seeds, Voronoi cells around 15 more
#' # regularly spaced seeds (Gibbs process), and random rectangles of 5 to 20
#' # cells a side; drawn in shuffled colours, since neighbouring region numbers
#' # would otherwise get similar ones
#' ts <- m |>
#'   drw_tessellation(type = "voronoi", n = 20, seed = 1, name = "voronoi") |>
#'   drw_tessellation(type = "gibbs", n = 15, interaction = 0.12, seed = 1,
#'                    name = "gibbs") |>
#'   drw_tessellation(type = "rectangle", size = c(5, 20), seed = 1,
#'                    name = "rectangle")
#' set.seed(1)
#' shuffled <- function(layer) {
#'   sample(hcl.colors(max(msk_pull(ts, layer)), "Dark 3"))
#' }
#' msk_vis(ts,
#'         .layer("voronoi", colours = shuffled("voronoi"), legend = FALSE),
#'         .layer("gibbs", colours = shuffled("gibbs"), legend = FALSE),
#'         .layer("rectangle", colours = shuffled("rectangle"), legend = FALSE))
#' @importFrom checkmate assertClass assertIntegerish assertChoice assertNumber
#'   assertCharacter
#' @export

drw_tessellation <- function(obj, type = "voronoi", n = 20L,
                              interaction = 0.1, size = c(5L, 20L),
                              name = "values", seed = NULL) {

  step <- .step()

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type, choices = c("voronoi", "rectangle", "gibbs"))
  assertIntegerish(x = n, len = 1, lower = 2)
  assertNumber(x = interaction, lower = 0, upper = 1)
  assertIntegerish(x = size, len = 2, lower = 1, any.missing = FALSE,
                   sorted = TRUE)
  assertCharacter(x = name, len = 1)
  assertIntegerish(x = seed, len = 1, null.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  extent <- obj@extent
  ncols <- obj@dims[1]
  nrows <- obj@dims[2]
  n <- as.integer(n)

  res_x <- (extent[2] - extent[1]) / ncols
  res_y <- (extent[4] - extent[3]) / nrows

  vals <- switch(type,

    voronoi = {
      sx <- stats::runif(n, extent[1], extent[2])
      sy <- stats::runif(n, extent[3], extent[4])
      sv <- seq_len(n)
      tesselateCpp(xmin = extent[1], ymax = extent[4],
                   res_x = res_x, res_y = res_y,
                   ncol = ncols, nrow = nrows,
                   seed_x = sx, seed_y = sy,
                   seed_val = as.integer(sv))
    },

    rectangle = {
      # dead leaves: drop rectangles at random until every cell is covered,
      # each on top of the earlier ones. A rectangle may hang over the edge,
      # so edge cells are as likely to be covered as inner ones
      grid <- matrix(0L, nrow = nrows, ncol = ncols)
      lens <- seq(size[1], size[2])
      drop <- 0L
      while (any(grid == 0L)) {
        drop <- drop + 1L
        w <- lens[sample.int(length(lens), 1)]
        h <- lens[sample.int(length(lens), 1)]
        x0 <- sample(seq(2L - w, ncols), 1)
        y0 <- sample(seq(2L - h, nrows), 1)
        cols <- max(1L, x0):min(ncols, x0 + w - 1L)
        rows <- max(1L, y0):min(nrows, y0 + h - 1L)
        grid[rows, cols] <- drop
      }
      # number the visible rectangles 1, 2, ... in row-major order
      v <- as.vector(t(grid))
      match(v, unique(v))
    },

    gibbs = {
      diag_len <- sqrt((extent[2] - extent[1])^2 + (extent[4] - extent[3])^2)
      r_inhib <- interaction * diag_len

      sx <- numeric(0)
      sy <- numeric(0)
      max_attempts <- n * 100L
      attempts <- 0L

      while (length(sx) < n && attempts < max_attempts) {
        cx <- stats::runif(1, extent[1], extent[2])
        cy <- stats::runif(1, extent[3], extent[4])
        attempts <- attempts + 1L

        if (length(sx) == 0) {
          sx <- cx
          sy <- cy
        } else {
          dists <- sqrt((sx - cx)^2 + (sy - cy)^2)
          if (all(dists >= r_inhib)) {
            sx <- c(sx, cx)
            sy <- c(sy, cy)
          }
        }
      }

      if (length(sx) < n) {
        warning("only ", length(sx), " of ", n, " seeds fit at this ",
                "'interaction'; reduce it or 'n' for all of them.",
                call. = FALSE)
      }
      if (length(sx) < 2L) {
        while (length(sx) < 2L) {
          sx <- c(sx, stats::runif(1, extent[1], extent[2]))
          sy <- c(sy, stats::runif(1, extent[3], extent[4]))
        }
      }

      sv <- seq_along(sx)
      tesselateCpp(xmin = extent[1], ymax = extent[4],
                   res_x = res_x, res_y = res_y,
                   ncol = ncols, nrow = nrows,
                   seed_x = sx, seed_y = sy,
                   seed_val = as.integer(sv))
    }
  )


  .update_mosaik(obj, values = as.numeric(vals), step = step)
}
