#' mosaik: Categorical Raster Algebra for Landscape Analysis
#'
#' Landscape-metric software has largely been a matter of fixed metric lists:
#' you get the metrics the author implemented, and a metric that is not on the
#' list is a feature request. mosaik takes the other route. It exposes the
#' \strong{measurement primitives and the means of composing them}, so that a
#' new metric is an expression rather than a new function.
#'
#' A landscape is a grid of stacked layers — the \code{\link{mosaik}} class.
#' Layers may be continuous fields (biomass, water, radiation), categorical
#' biophysical states (cover, soil, vegetation stage), or categorical
#' institutional states (designation, tenure). All are held in the same data
#' model and every operation treats them identically, whether they were measured
#' or synthesised.
#'
#' Four things distinguish it:
#'
#' \describe{
#'   \item{Categorical-first}{Cell values \emph{are} group IDs, and
#'     \code{@categories} and \code{@patches} are first-class slots rather than
#'     levels bolted onto a numeric raster.}
#'   \item{Native provenance}{\code{@provenance} is a slot that survives every
#'     operation, not a log written alongside the result.}
#'   \item{Composability}{Mosaic in, mosaic out, so operations chain with the
#'     pipe; and \code{\link{msr}} derives a metric from an equation over the
#'     primitives in \code{metric.scale} notation.}
#'   \item{Five primitives}{Adjacency, area, perimeter, distance and
#'     dissimilarity, from which nearly every published landscape metric can be
#'     composed.}
#' }
#'
#' The MSPA vignette is the demonstration: a published segmentation algorithm
#' reproduced at full fidelity out of the primitives alone.
#'
#' There are four function families:
#'
#' \describe{
#'   \item{\code{syn_*}}{synthesise an abstract field — noise, texture,
#'     gradients, point patterns, clusters, tessellations
#'     (\code{\link{syn_noise}}, \code{\link{syn_texture}},
#'     \code{\link{syn_gradient}}, \code{\link{syn_pattern}},
#'     \code{\link{syn_cluster}}, \code{\link{syn_tessellation}}). These make
#'     fields to analyse, not landscapes; for landscapes with a known ground
#'     truth see the \pkg{mundus} package.}
#'   \item{\code{mdf_*}}{modify layers with generic operators that claim no new
#'     semantics (e.g. \code{\link{mdf_morph}}, \code{\link{mdf_distance}},
#'     \code{\link{mdf_fill}}, \code{\link{mdf_scale}}). A recorded sequence of
#'     them is a recipe, replayable on any raster with \code{\link{mdf}}.}
#'   \item{\code{msr_*} and \code{\link{msr}}}{measure — \code{msr_*} for the
#'     primitives, \code{\link{msr}} for a metric derived from an equation over
#'     them.}
#'   \item{\code{msk_*}}{utilities and accessors (\code{\link{msk_vis}},
#'     \code{\link{msk_extent}}, \code{\link{msk_terra}}, ...).}
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
