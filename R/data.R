#' Example landscape
#'
#' A 60 x 56 artificial landscape with two layers: a categorical land-cover map
#' and a continuous intensity surface. The data originate from a hand-crafted
#' landscape mosaic used during early package development.
#'
#' @format A \code{mosaik} object with extent \code{c(0, 60, 0, 56)},
#'   resolution 1, and two layers:
#' \describe{
#'   \item{cover}{Integer land-cover class (9 classes: 1, 11, 21, 24, 27, 31,
#'     41, 44, 47).}
#'   \item{intensity}{Continuous values in 0--100.}
#' }
#' @source Manually constructed for illustration purposes.
"landscape"

#' Binary pattern for morphological analysis
#'
#' The 38 x 31 binary pattern of Fig. 1 in Soille & Vogt (2009), the example
#' that paper uses to introduce morphological spatial pattern analysis.
#' \code{vignette("mspa")} builds the full seven-class segmentation of it from
#' ordinary operators.
#'
#' @format A \code{mosaik} object with extent \code{c(0, 38, 0, 31)},
#'   resolution 1, and one layer:
#' \describe{
#'   \item{pattern}{1 foreground (habitat), 0 background (matrix). 710 of the
#'     1178 cells are foreground.}
#' }
#' @references Soille, P. & Vogt, P. (2009). Morphological segmentation of
#'   binary patterns. \emph{Pattern Recognition Letters} 30(4), 456-459.
#' @examples
#' msk_dims(mspa)
#' sum(msk_pull(mspa, "pattern"))
"mspa"
