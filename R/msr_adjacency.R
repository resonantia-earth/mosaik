#' Adjacency of classes
#'
#' Count how often the cells of each class of a layer border cells of the same
#' or another class.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param connect [`integerish(1)`][integer]\cr which cells border each other:
#'   \code{4} (default) the orthogonal neighbours, \code{8} also the diagonal
#'   ones.
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with the class by class matrix \code{adjacency}
#'   added to the class table of \code{layer} (see
#'   \code{\link{msk_categories}}).
#' @details \code{adjacency} counts the pairs of bordering cells. Each pair is
#'   counted from both sides, as in FRAGSTATS: two forest cells side by side
#'   add 2 to forest-forest, a forest cell beside a meadow cell adds 1 to
#'   forest-meadow and 1 to meadow-forest. The matrix is therefore symmetric,
#'   and its diagonal holds how often a class borders itself. In
#'   \code{\link{msr}}, with a class in focus, \code{adjacency.self} is that
#'   class's own adjacencies and \code{adjacency.self + sum(adjacency.others)}
#'   all adjacencies of the class.
#' @examples
#' # how often each land cover class borders each other class
#' m <- msr_adjacency(landscape, layer = "cover")
#' msk_categories(m, layer = "cover")$adjacency
#'
#' # where the forest patches touch each other, with diagonal neighbours
#' f <- mdf_filter(landscape, cover == 47, add = "forest") |>
#'   mdf_componentise(connectivity = 4L, layer = "forest", add = "patch") |>
#'   msr_adjacency(connect = 8, layer = "patch")
#' msk_categories(f, layer = "patch")$adjacency
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_adjacency <- function(obj = NULL, connect = 4, layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = connect, choices = c(4, 8))
  assertCharacter(x = layer, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims

  uVals <- sort(unique(vals[!is.na(vals)]))

  # every pair of bordering cells counted from both sides
  adjacency <- countCellAdjacenciesCpp(vals = vals, nrow = dims[2],
                                       ncol = dims[1], doublecount = TRUE,
                                       eightconn = connect == 8)

  obj <- .store_class(obj, layer, "adjacency", as.integer(uVals), adjacency)

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
