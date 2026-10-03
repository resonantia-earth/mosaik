#' Make a textured surface
#'
#' Generate multi-scale or statistically structured continuous surfaces. These
#' produce spatially autocorrelated patterns with controlled roughness.
#'
#' @param obj [`mosaik`]\cr the mosaik whose grid the texture is generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [`character(1)`][character]\cr the texture type:
#'   \code{"diamondSquare"} (default), \code{"fbm"} (fractional Brownian
#'   motion), \code{"billow"} (abs-value fBm, cloud-like lumps), or
#'   \code{"ridged"} (inverted ridged fBm).
#' @param base [`character(1)`][character]\cr base noise function for
#'   fBm/billow/ridged octave stacking: \code{"perlin"} (default) or
#'   \code{"simplex"}. Ignored for \code{"diamondSquare"}.
#' @param hurst [`numeric(1)`][numeric]\cr the Hurst exponent controlling
#'   roughness (0 = very rough, 1 = very smooth). Default 0.7.
#' @param octaves [`integerish(1)`][integer]\cr number of octaves for
#'   fBm/billow/ridged. Default 6. More octaves add finer detail.
#' @param lacunarity [`numeric(1)`][numeric]\cr frequency multiplier per octave.
#'   Default 2.0.
#' @param frequency [`numeric(1)`][numeric]\cr base frequency for the first
#'   octave. Default 4.
#' @param startDev [`numeric(1)`][numeric]\cr initial standard deviation for
#'   diamond-square. Default 1.
#' @param name [`character(1)`][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param seed [`integerish(1)`][integer]\cr random seed for reproducibility.
#' @return A \code{\link{mosaik}} with texture values scaled to \[0, 1\].
#' @references Fournier A, Fussell D, Carpenter L. Computer rendering of
#'   stochastic models. Communications of the ACM. 1982;25:371-384.
#'
#'   Mandelbrot BB, Van Ness JW. Fractional Brownian motions, fractional noises
#'   and applications. SIAM Review. 1968;10(4):422-437.
#'
#'   Perlin K. Improving noise. ACM Transactions on Graphics (SIGGRAPH 2002).
#'   2002;21(3):681-682.
#'
#'   Ebert DS, Musgrave FK, Peachey D, Perlin K, Worley S. Texturing and
#'   Modeling: A Procedural Approach. 3rd ed. Morgan Kaufmann; 2003.
#' @family generator functions
#' @examples
#' m <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # diamond-square heightmap, rough and smooth
#' t <- m |>
#'   drw_texture(type = "diamondSquare", hurst = 0.3, seed = 1, name = "rough") |>
#'   drw_texture(type = "diamondSquare", hurst = 0.9, seed = 1, name = "smooth")
#' msk_vis(t, .layer("rough"), .layer("smooth"))
#'
#' # fBm, billow (cloud-like lumps) and ridged (mountain ridges) on Perlin noise
#' t <- m |>
#'   drw_texture(type = "fbm", seed = 1, name = "fbm") |>
#'   drw_texture(type = "billow", seed = 1, name = "billow") |>
#'   drw_texture(type = "ridged", seed = 1, name = "ridged")
#' msk_vis(t, .layer("fbm"), .layer("billow"), .layer("ridged"))
#'
#' # the same on Simplex noise
#' t <- m |>
#'   drw_texture(type = "fbm", base = "simplex", seed = 1, name = "fbm") |>
#'   drw_texture(type = "billow", base = "simplex", seed = 1, name = "billow") |>
#'   drw_texture(type = "ridged", base = "simplex", seed = 1, name = "ridged")
#' msk_vis(t, .layer("fbm"), .layer("billow"), .layer("ridged"))
#' @importFrom checkmate assertClass assertIntegerish assertChoice assertNumber
#'   assertCharacter
#' @export

drw_texture <- function(obj, type = "diamondSquare", base = "perlin",
                         hurst = 0.7, octaves = 6L, lacunarity = 2.0,
                         frequency = 4, startDev = 1, name = "values",
                         seed = NULL) {

  step <- .step()

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type, choices = c("diamondSquare", "fbm", "billow", "ridged"))
  assertChoice(x = base, choices = c("perlin", "simplex"))
  assertNumber(x = hurst, lower = 0, upper = 1)
  assertIntegerish(x = octaves, len = 1, lower = 1, upper = 16)
  assertNumber(x = lacunarity, lower = 1)
  assertNumber(x = frequency, lower = 0.01)
  assertNumber(x = startDev, lower = 0)
  assertCharacter(x = name, len = 1)
  assertIntegerish(x = seed, len = 1, null.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  ncols <- obj@dims[1]
  nrows <- obj@dims[2]
  n_cells <- ncols * nrows

  vals <- switch(type,

    diamondSquare = {
      k <- ceiling(log2(max(ncols, nrows)))
      size <- 2L^k + 1L
      mat <- matrix(0, nrow = size, ncol = size)

      mat[1, 1]       <- stats::rnorm(1, 0, startDev)
      mat[1, size]    <- stats::rnorm(1, 0, startDev)
      mat[size, 1]    <- stats::rnorm(1, 0, startDev)
      mat[size, size] <- stats::rnorm(1, 0, startDev)

      stepSizes <- 2L^(k:1)
      roughness <- if (is.null(hurst)) 0.5 else hurst

      result <- diamondSquareCpp(mat = mat, stepSize = stepSizes,
                                 roughness = roughness, startDev = startDev)

      cropped <- result[seq_len(nrows), seq_len(ncols)]
      as.vector(t(cropped))
    },

    fbm = {
      .fbm_octaves(ncols, nrows, base, frequency, octaves, lacunarity,
                    hurst, transform = "none")
    },

    billow = {
      .fbm_octaves(ncols, nrows, base, frequency, octaves, lacunarity,
                    hurst, transform = "billow")
    },

    ridged = {
      .fbm_octaves(ncols, nrows, base, frequency, octaves, lacunarity,
                    hurst, transform = "ridged")
    }
  )

  # scale to [0, 1]
  rng <- range(vals, na.rm = TRUE)
  if (rng[2] > rng[1]) {
    vals <- (vals - rng[1]) / (rng[2] - rng[1])
  } else {
    vals <- rep(0.5, n_cells)
  }


  .update_mosaik(obj, values = vals, step = step)
}


#' Stack octaves of base noise into fractional Brownian motion
#'
#' Sums \code{octaves} layers of Perlin/Simplex noise, each at a higher
#' frequency (multiplied by \code{lacunarity}) and lower amplitude than the
#' last. The amplitude gain per octave is \code{2^(-hurst)}: a higher Hurst
#' exponent decays amplitude faster, yielding a smoother surface. Each octave
#' is offset by a fresh random shift so octaves are decorrelated; the caller
#' is responsible for \code{set.seed()}.
#'
#' @param ncols,nrows grid dimensions.
#' @param base \code{"perlin"} or \code{"simplex"}.
#' @param frequency base frequency for the first octave.
#' @param octaves number of octaves to stack.
#' @param lacunarity frequency multiplier per octave.
#' @param hurst Hurst exponent in \[0, 1\]; amplitude gain is \code{2^(-hurst)}.
#' @param transform \code{"none"} (plain fBm), \code{"billow"}
#'   (\code{abs}, cloud-like lumps), or \code{"ridged"} (\code{1 - abs},
#'   mountain ridges) applied to each octave before summation.
#' @return A flat numeric vector of length \code{ncols * nrows} (row-major),
#'   unscaled — \code{drw_texture} rescales to \[0, 1\].
#' @noRd

.fbm_octaves <- function(ncols, nrows, base, frequency, octaves, lacunarity,
                         hurst, transform = "none") {

  noise_fn <- if (base == "simplex") simplexCpp else perlinCpp
  gain <- 2^(-hurst)
  n_cells <- ncols * nrows

  total <- numeric(n_cells)
  amp <- 1
  freq <- frequency

  for (o in seq_len(octaves)) {
    offset_x <- stats::runif(1, 0, 1000)
    offset_y <- stats::runif(1, 0, 1000)
    layer <- noise_fn(ncol = ncols, nrow = nrows, freq = freq,
                      offset_x = offset_x, offset_y = offset_y)

    layer <- switch(transform,
      none   = layer,
      billow = abs(layer),
      ridged = 1 - abs(layer)
    )

    total <- total + amp * layer
    amp <- amp * gain
    freq <- freq * lacunarity
  }

  total
}
