#' mosaik: Categorical Raster Algebra for Landscape Analysis
#'
#' Landscape tools usually come as fixed lists: a catalogue of metrics, a
#' programme for one pattern analysis, a generator for one kind of neutral
#' landscape. mosaik offers the \strong{parts those tools are made of and the
#' means of combining them} instead: operations that synthesise fields, modify
#' layers and measure them. A new metric, a new morphological analysis or a new
#' synthetic test field is then a combination of existing operations, not a new
#' function.
#'
#' Combining operations only works if they share what they know: an equation
#' over the perimeter and the area of a patch needs both to refer to the same
#' patch, a step that measures, modifies and measures again needs the earlier
#' results to still be there, and a reported number needs to be traceable
#' through every step that made it. The \code{\link{mosaik}} class therefore
#' holds the layers, the patches, the classes, the results and the history of
#' operations, and every operation keeps them consistent.
#'
#' A mosaik can be a grid of stacked layers, categorical (land cover, soil,
#' tenure) or continuous (biomass, intensity). Cell values of a categorical
#' layer are the class identifiers themselves. Results are written into the
#' object at the level they describe: patch-level in \code{@patches},
#' class-level in \code{@categories}, landscape-level in \code{@global}.
#'
#' The primitives are area, number, perimeter, adjacency, dissimilarity and
#' cost (distance is the cost measured in metres), each at patch, class or
#' landscape level. The MSPA vignette shows how far they reach: a published
#' segmentation algorithm rebuilt from them, matching the original cell for
#' cell.
#'
#' There are four function families:
#'
#' \describe{
#'   \item{\code{syn_*}}{synthesise an abstract field to analyse: noise,
#'     texture, gradients, patterns, clusters, tessellations
#'     (\code{\link{syn_noise}}, \code{\link{syn_texture}},
#'     \code{\link{syn_gradient}}, \code{\link{syn_pattern}},
#'     \code{\link{syn_cluster}}, \code{\link{syn_tessellation}}).}
#'   \item{\code{mdf_*}}{modify layers with generic operators that do not know
#'     what a layer represents (e.g. \code{\link{mdf_morph}},
#'     \code{\link{mdf_distance}}, \code{\link{mdf_fill}},
#'     \code{\link{mdf_scale}}). A recorded sequence of them is a recipe,
#'     replayable on any raster with \code{\link{mdf}}.}
#'   \item{\code{msr_*} and \code{\link{msr}}}{measure: \code{msr_*} for the
#'     primitives, \code{\link{msr}} for a metric derived from an equation over
#'     them in \code{metric.scale} notation.}
#'   \item{\code{msk_*}}{accessors and utilities (\code{\link{msk_vis}},
#'     \code{\link{msk_categories}}, \code{\link{msk_terra}}, ...).}
#' }
#'
#' @author \strong{Maintainer, Author}: Steffen Ehrmann
#'   \email{steffen.ehrmann@posteo.de}
#'
#' @seealso \itemize{ \item Github project:
#'   \href{https://github.com/resonantia-earth/mosaik}{https://github.com/resonantia-earth/mosaik}
#'   \item Report bugs:
#'   \href{https://github.com/resonantia-earth/mosaik/issues}{https://github.com/resonantia-earth/mosaik/issues}
#'   }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib mosaik, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
