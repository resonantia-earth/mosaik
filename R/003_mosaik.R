#' Create a landscape mosaik
#'
#' \code{mosaik} creates a gridded landscape mosaik, either from scratch or by
#' importing a \code{SpatRaster} (requires the \pkg{terra} package).
#'
#' @param extent [numeric(4)][numeric]\cr the extent as \code{c(xmin, xmax,
#'   ymin, ymax)}. Ignored when \code{rast} is supplied.
#' @param res [numeric(1) | numeric(2)][numeric]\cr the cell size (resolution).
#'   Recycled if length 1. Ignored when \code{rast} is supplied.
#' @param crs [character(1)][character]\cr coordinate reference system (proj4
#'   string), or \code{NA_character_} for cartesian (the default). Ignored when
#'   \code{rast} is supplied (CRS is read from the raster).
#' @param vals a matrix, a named list of matrices or vectors, or \code{NULL}
#'   (default) to create an empty grid filled with \code{NA}. Ignored when
#'   \code{rast} is supplied.
#' @param rast a \code{SpatRaster} object to import. When supplied,
#'   \code{extent}, \code{res}, \code{crs}, and \code{vals} are ignored.
#' @param group [logical(1)][logical]\cr if \code{TRUE} and the raster has no
#'   categories, unique values are stored as categories. Only used when
#'   \code{rast} is supplied.
#' @return An object of class \code{\link{mosaik-class}}.
#' @examples
#' # create an empty 100x100 mosaik with unit cells
#' msk <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # ... with finer resolution and a geographic CRS
#' msk <- mosaik(extent = c(5, 15, 47, 55), res = 0.01,
#'               crs = "+proj=longlat +datum=WGS84")
#'
#' # ... from a matrix
#' msk <- mosaik(extent = c(0, 20, 0, 20), res = 1,
#'               vals = matrix(runif(400), nrow = 20, ncol = 20))
#'
#' # ... from named layers
#' msk <-  mosaik(extent = c(0, 10, 0, 10), res = 1,
#'                vals = list(cover = sample(1:5, 100, replace = TRUE),
#'                            ndvi = runif(100)))
#'
#' # ... from a SpatRaster
#' \dontrun{
#' library(terra)
#' r <- rast(ncols = 10, nrows = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
#' values(r) <- sample(1:5, 100, replace = TRUE)
#' msk <- mosaik(rast = r)
#' }
#' @importFrom checkmate assertNumeric assertCharacter assertLogical
#' @export

mosaik <- function(extent = NULL, res = NULL, crs = NA_character_,
                   vals = NULL, rast = NULL, group = FALSE) {

  assertLogical(x = group, len = 1)

  # Each input path fills these slots, differently; the object is constructed
  # once, at the end, so new_mosaik()'s contract is honoured in exactly one place.
  categories <- list()

  if (!is.null(rast)) {

    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' must be installed to import a SpatRaster.",
           call. = FALSE)
    }

    extent <- as.vector(terra::ext(rast))
    dims <- c(terra::ncol(rast), terra::nrow(rast))
    crs <- as.character(terra::crs(rast, proj = TRUE))
    if (is.null(crs) || crs == "") crs <- NA_character_

    layers <- list()
    provenance <- list()
    for (i in seq_along(names(rast))) {
      nm <- names(rast)[i]
      rawVal <- terra::values(rast[[i]])[, 1]
      layers[[nm]] <- rawVal
      provenance <- c(provenance,
                      list(paste0("layer '", nm, "' imported from SpatRaster.")))

      cats <- terra::cats(rast[[i]])
      if (is.data.frame(cats[[1]]) && nrow(cats[[1]]) > 0) {
        catDF <- cats[[1]]
        gid <- catDF[[1]]
        val <- if (ncol(catDF) > 1) catDF[[2]] else as.character(gid)
        categories[[nm]] <- list(gid = gid, val = val)
      } else if (group) {
        uvals <- sort(unique(rawVal[!is.na(rawVal)]))
        categories[[nm]] <- list(gid = uvals, val = as.character(uvals))
      }
    }

  } else {

    if (is.null(extent)) stop("'extent' is required when 'rast' is not supplied.")
    if (is.null(res)) stop("'res' is required when 'rast' is not supplied.")

    assertNumeric(x = extent, len = 4, any.missing = FALSE)
    assertNumeric(x = res, min.len = 1, max.len = 2, lower = 0,
                  any.missing = FALSE)
    assertCharacter(x = crs, len = 1)

    if (length(res) == 1) res <- rep(res, 2)
    dims <- c(as.integer(round((extent[2] - extent[1]) / res[1])),
              as.integer(round((extent[4] - extent[3]) / res[2])))

    if (is.null(vals)) {
      layers <- list()
    } else if (is.matrix(vals)) {
      layers <- list(values = as.vector(t(vals)))
    } else if (is.list(vals)) {
      layers <- lapply(vals, function(v) {
        if (is.matrix(v)) as.vector(t(v)) else v
      })
    } else {
      stop("'vals' must be NULL, a matrix, or a named list.")
    }

    provenance <- list("created via mosaik()")
  }

  new_mosaik(
    extent     = extent,
    dims       = dims,
    layers     = layers,
    categories = categories,
    crs        = crs,
    provenance = provenance
  )
}
