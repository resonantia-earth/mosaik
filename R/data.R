#' Example landscape
#'
#' A 60 x 56 artificial landscape with a land-cover map and a map of canopy
#' height. Its structures are placed deliberately, so that the analyses in
#' the vignettes find something known: forest patches of very different sizes,
#' an enclosed glade and an old windthrow, a fresh clear-cut cut in from a patch's edge, a corridor
#' between two cores and one interrupted by a river, a dead-end branch and a
#' loop, a road that runs beside the river and then cuts through a forest, a
#' decaying bocage of hedgerows, sharp and gradual forest edges, and patches
#' bordering few or many other classes.
#'
#' @format A \code{mosaik} object with extent \code{c(0, 60, 0, 56)},
#'   resolution 1, and two layers:
#' \describe{
#'   \item{cover}{Land-cover class, labelled in the categories of the layer:
#'     1 river, 11 arable land, 21 intensive grassland, 24 extensive
#'     grassland, 27 fallow, shrub and clear-cuts, 31 settlement, 35 road,
#'     41 wetland, 44 orchard, 47 forest and hedgerows.}
#'   \item{canopy}{Canopy height in metres. Old forest stands about 25 m
#'     tall and falls towards its edges, except along the road and the fresh
#'     clear-cuts, where the edge is sharp; the forest in the south-west is
#'     younger and lower. Outside the forest it shows hedgerows, gardens,
#'     orchards and single trees in the fields, which the land cover does not.}
#' }
#' @source Built by \code{data-raw/landscape.R} in the package sources.
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
