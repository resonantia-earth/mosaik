# TODO: this seems least developed. Not sure those patterns make much sense.
# Please pull the nlmpy/NLMR paper and check the algos there.

#' Make a tessellation surface
#'
#' Partition the grid into discrete regions. \code{"voronoi"} places random seed
#' points and assigns each cell to the nearest seed. \code{"rectangle"}
#' generates random rectangular blocks. \code{"gibbs"} uses a Strauss
#' interaction process for more regular seed spacing before Voronoi assignment.
#'
#' @param obj [mosaik]\cr the mosaik whose grid the tessellation is generated
#'   on. Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [character(1)][character]\cr tessellation type: \code{"voronoi"}
#'   (default), \code{"rectangle"}, or \code{"gibbs"}.
#' @param n [integerish(1)][integer]\cr number of regions (seed points for
#'   voronoi/gibbs, rectangular blocks for rectangle). Default 20.
#' @param interaction [numeric(1)][numeric]\cr for \code{"gibbs"} only: the
#'   inhibition radius as a fraction of the grid diagonal. Points closer than
#'   this are penalised. Default 0.1.
#' @param name [character(1)][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param seed [integerish(1)][integer]\cr random seed for reproducibility.
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
#' # Voronoi tessellation with 20 random seeds
#' syn_tessellation(m, type = "voronoi", n = 20, seed = 1)
#'
#' # random rectangular blocks
#' syn_tessellation(m, type = "rectangle", n = 12, seed = 1)
#'
#' # Gibbs process, more regular spacing
#' syn_tessellation(m, type = "gibbs", n = 15, interaction = 0.12, seed = 1)
#' @importFrom checkmate assertClass assertIntegerish assertChoice assertNumber
#'   assertCharacter
#' @export

syn_tessellation <- function(obj, type = "voronoi", n = 20L,
                              interaction = 0.1, name = "values",
                              seed = NULL) {

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type, choices = c("voronoi", "rectangle", "gibbs"))
  assertIntegerish(x = n, len = 1, lower = 2)
  assertNumber(x = interaction, lower = 0, upper = 1)
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
      out <- rep(0L, ncols * nrows)
      col_idx <- rep(seq_len(ncols), times = nrows) - 1L
      row_idx <- rep(seq_len(nrows), each = ncols) - 1L

      n_splits <- max(2L, as.integer(ceiling(sqrt(n))))
      x_breaks <- sort(c(0, sample.int(ncols - 1L, n_splits - 1L), ncols))
      y_breaks <- sort(c(0, sample.int(nrows - 1L, n_splits - 1L), nrows))

      label <- 1L
      for (i in seq_len(length(x_breaks) - 1)) {
        for (j in seq_len(length(y_breaks) - 1)) {
          mask <- col_idx >= x_breaks[i] & col_idx < x_breaks[i + 1] &
                  row_idx >= y_breaks[j] & row_idx < y_breaks[j + 1]
          out[mask] <- label
          label <- label + 1L
        }
      }
      out
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

      if (length(sx) < 2L) {
        warning("Gibbs process could only place ", length(sx),
                " seeds. Try reducing 'interaction' or 'n'.")
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

  prov <- msk_prov("syn_tessellation",
                     list(type = type, n = n, interaction = interaction,
                          name = name, seed = seed))

  msk_set(obj, name, as.numeric(vals), prov)
}
