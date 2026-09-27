# TODO (design notes, not code):
# - probably better to call this mdf_rasterize/mdf_rasterise
# - what does this function provide over fasterize? Test whether it is similarly
#   fast; if not, revise to make it so, and note that this also records
#   provenance on top.

#' Burn polygon zones into a mosaik layer
#'
#' Rasterise one or more polygons into a mosaik layer by testing cell centroids
#' against polygon boundaries. Each polygon assigns its \code{value} to every
#' cell whose centroid falls inside. Later polygons overwrite earlier ones where
#' they overlap.
#'
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param geom either a single closed polygon as a two-column numeric matrix
#'   (x, y; first row == last row), or a \code{list} of such matrices for
#'   multiple polygons.
#' @param value [integerish(1)][integer]\cr the value(s) to assign to cells inside the
#'   polygon(s). A single value is recycled for all polygons; a vector must
#'   match the number of polygons.
#' @param background [integerish(1)][integer]\cr value for cells outside all
#'   polygons. Defaults to \code{NA}. Set to \code{NULL} to leave existing layer
#'   values untouched outside polygons.
#' @param layer [character(1)][character]\cr the layer to write into. Defaults
#'   to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj} with polygon zones
#'   burned into the specified layer.
#' @examples
#' # single polygon zone
#' poly <- matrix(c(10,10, 50,10, 50,40, 10,40, 10,10),
#'                ncol = 2, byrow = TRUE)
#' mdf_zonify(landscape, geom = poly, value = 99L, layer = "cover")
#'
#' # multiple polygon zones
#' p1 <- matrix(c(5,5, 25,5, 25,25, 5,25, 5,5), ncol = 2, byrow = TRUE)
#' p2 <- matrix(c(35,30, 55,30, 55,50, 35,50, 35,30), ncol = 2, byrow = TRUE)
#' mdf_zonify(landscape, geom = list(p1, p2), value = c(88L, 99L),
#'            layer = "cover")
#'
#' # overlay on existing values (background = NULL preserves original)
#' mdf_zonify(landscape, geom = poly, value = 99L, background = NULL,
#'            layer = "cover")
#' @family operators to select a subset of cells
#' @importFrom checkmate assertClass assertIntegerish assertCharacter assertMatrix assertList
#' @export

mdf_zonify <- function(obj = NULL,
                       geom,
                       value,
                       background = NA,
                       layer = NULL,
                       add = NULL) {

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # normalise geom to a list of matrices
  if (is.matrix(geom)) {
    geom <- list(geom)
  }
  assertList(x = geom, types = "matrix", min.len = 1)

  for (i in seq_along(geom)) {
    g <- geom[[i]]
    if (ncol(g) != 2) stop("Each polygon matrix must have 2 columns (x, y).")
    if (nrow(g) < 4) stop("Each polygon must have at least 4 rows (3 vertices + closing).")
    if (g[1, 1] != g[nrow(g), 1] || g[1, 2] != g[nrow(g), 2]) {
      stop("Polygon ", i, " is not closed (first row must equal last row).")
    }
  }

  # recycle or validate value
  assertIntegerish(x = value, min.len = 1)
  if (length(value) == 1) {
    value <- rep(as.integer(value), length(geom))
  } else if (length(value) != length(geom)) {
    stop("'value' must be length 1 or match the number of polygons (", length(geom), ").")
  }

  # pull data ----
  if (is.null(layer)) layer <- names(obj@layers)[1]
  dims <- obj@dims
  ext  <- obj@extent
  res  <- msk_res(obj)

  # build cell centroid coordinates
  xs <- seq(ext[1] + res[1] / 2, by = res[1], length.out = dims[1])
  ys <- seq(ext[4] - res[2] / 2, by = -res[2], length.out = dims[2])
  coords_x <- rep(xs, times = dims[2])
  coords_y <- rep(ys, each = dims[1])
  pts <- cbind(coords_x, coords_y)

  # initialise output
  if (is.null(background)) {
    temp <- msk_pull(obj, layer)
  } else {
    temp <- rep(background, prod(dims))
  }

  # burn each polygon ----
  for (i in seq_along(geom)) {
    inside <- pointInPolyCpp(vert = pts, geom = geom[[i]], invert = FALSE)
    temp[inside >= 1L] <- value[i]
  }

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_zonify", list(n_polygons = length(geom),
                                        values = value,
                                        background = background,
                                        layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
