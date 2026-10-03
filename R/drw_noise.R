# TODO: this looks right, but is there any way to validate these types?

#' Make a noise surface
#'
#' Generate spatially uncorrelated or single-scale noise surfaces. These are
#' building blocks for more complex textures (see \code{\link{drw_texture}}).
#'
#' @param obj [`mosaik`]\cr the mosaik whose grid the noise is generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [`character(1)`][character]\cr the noise type:
#'   \code{"white"} (default), \code{"perlin"}, or \code{"simplex"}.
#' @param frequency [`numeric(1)`][numeric]\cr spatial frequency for structured
#'   noise types (Perlin, Simplex). Higher values produce finer detail.
#'   Default 4.
#' @param name [`character(1)`][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param seed [`integerish(1)`][integer]\cr random seed for reproducibility.
#' @return A \code{\link{mosaik}} with noise values scaled to \[0, 1\].
#' @references
#'   Perlin K. Improving noise. ACM Transactions on Graphics (SIGGRAPH 2002).
#'   2002;21(3):681-682.
#'
#'   Perlin K. Noise hardware. In: Olano M (ed.), Real-Time Shading SIGGRAPH
#'   Course Notes. 2001. (Simplex noise)
#'
#'   Gustavson S. Simplex noise demystified. Linkoping University, Sweden. 2005.
#' @family generator functions
#' @examples
#' m <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # uniform random noise, Perlin noise and Simplex noise
#' n <- m |>
#'   drw_noise(type = "white", seed = 42, name = "white") |>
#'   drw_noise(type = "perlin", frequency = 3, seed = 1, name = "perlin") |>
#'   drw_noise(type = "simplex", frequency = 3, seed = 1, name = "simplex")
#' msk_vis(n, .layer("white"), .layer("perlin"), .layer("simplex"))
#'
#' # a higher frequency gives finer detail
#' p <- m |>
#'   drw_noise(type = "perlin", frequency = 2, seed = 1, name = "frequency_2") |>
#'   drw_noise(type = "perlin", frequency = 8, seed = 1, name = "frequency_8")
#' msk_vis(p, .layer("frequency_2"), .layer("frequency_8"))
#' @importFrom checkmate assertClass assertChoice assertNumber assertCharacter
#'   assertIntegerish
#' @export

drw_noise <- function(obj, type = "white", frequency = 4,
                       name = "values", seed = NULL) {

  step <- .step()

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type, choices = c("white", "perlin", "simplex"))
  assertNumber(x = frequency, lower = 0.01)
  assertCharacter(x = name, len = 1)
  assertIntegerish(x = seed, len = 1, null.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  ncols <- obj@dims[1]
  nrows <- obj@dims[2]
  n_cells <- ncols * nrows

  vals <- switch(type,

    white = {
      stats::runif(n_cells)
    },

    perlin = {
      offset_x <- stats::runif(1, 0, 1000)
      offset_y <- stats::runif(1, 0, 1000)
      raw <- perlinCpp(ncol = ncols, nrow = nrows, freq = frequency,
                       offset_x = offset_x, offset_y = offset_y)
      rng <- range(raw)
      if (rng[2] > rng[1]) (raw - rng[1]) / (rng[2] - rng[1]) else rep(0.5, n_cells)
    },

    simplex = {
      offset_x <- stats::runif(1, 0, 1000)
      offset_y <- stats::runif(1, 0, 1000)
      raw <- simplexCpp(ncol = ncols, nrow = nrows, freq = frequency,
                        offset_x = offset_x, offset_y = offset_y)
      rng <- range(raw)
      if (rng[2] > rng[1]) (raw - rng[1]) / (rng[2] - rng[1]) else rep(0.5, n_cells)
    }
  )


  .update_mosaik(obj, values = vals, step = step)
}
