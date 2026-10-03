#' Accessor functions for mosaik objects
#'
#' Simple functions to access properties of a \code{\link{mosaik}} object.
#'
#' @param x a \code{\link{mosaik}} object.
#' @param layer [`character(1)`][character]\cr for \code{msk_categories},
#'   \code{msk_patches} and \code{msk_global}: the layer whose results to
#'   return. Defaults to the first layer.
#' @return For \code{msk_categories}, \code{msk_patches} and
#'   \code{msk_global}: the class, patch or landscape-level results of
#'   \code{layer}, as a list, or \code{NULL} if none have been measured.
#' @name accessors
#' @examples
#' msk <- mosaik(extent = c(0, 10, 0, 10), res = 1,
#'               vals = list(cover = sample(1:5, 100, replace = TRUE)))
#' msk_extent(msk)
#' msk_dims(msk)
#' msk_res(msk)
#' msk_ncells(msk)
#' msk_crs(msk)
#' msk_names(msk)
#'
#' msk <- msr_area(msk, scale = "class", layer = "cover")
#' msk_categories(msk, layer = "cover")
NULL

#' @rdname accessors
#' @export
msk_extent <- function(x) x@extent

#' @rdname accessors
#' @export
msk_dims <- function(x) x@dims

#' @rdname accessors
#' @export
msk_res <- function(x) {
  c((x@extent[2] - x@extent[1]) / x@dims[1],
    (x@extent[4] - x@extent[3]) / x@dims[2])
}

#' @rdname accessors
#' @export
msk_ncells <- function(x) prod(x@dims)

#' @rdname accessors
#' @export
msk_crs <- function(x) x@crs

#' @rdname accessors
#' @export
msk_provenance <- function(x) x@provenance

#' @rdname accessors
#' @export
msk_names <- function(x) names(x@layers)

#' @rdname accessors
#' @export
msk_categories <- function(x, layer = NULL) .layer_entry(x, "categories", layer)

#' @rdname accessors
#' @export
msk_patches <- function(x, layer = NULL) .layer_entry(x, "patches", layer)

#' @rdname accessors
#' @export
msk_global <- function(x, layer = NULL) .layer_entry(x, "global", layer)

# the entry of one layer in a per-layer slot; NULL if nothing is stored for it
.layer_entry <- function(x, slot, layer) {
  if (is.null(layer)) layer <- names(x@layers)[1]
  if (!layer %in% names(x@layers)) {
    stop("layer '", layer, "' not found in 'x'.", call. = FALSE)
  }
  methods::slot(x, slot)[[layer]]
}
