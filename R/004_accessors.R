#' Accessor functions for mosaik objects
#'
#' Simple functions to access properties of a \code{\link{mosaik}} object.
#'
#' @param x a \code{\link{mosaik}} object.
#' @details What is measured on a layer is read with \code{\link{msk_table}}.
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
