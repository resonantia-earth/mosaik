#' Adjacency of cells
#'
#' Count how often the cells of each class border cells of the same or another
#' class, or how the patches of a layer touch each other.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param scale [`character(1)`][character]\cr \code{"class"} (default) counts
#'   the adjacencies between classes; \code{"patch"} those between the patches
#'   of \code{layer}, which must have been numbered with
#'   \code{\link{mdf_componentise}} first.
#' @param type [`character(1)`][character]\cr at class scale, which adjacencies
#'   to count: \code{"like"} (how often a class borders itself, a measure of
#'   clumping), \code{"paired"} (how often each class borders each other class)
#'   or \code{"pairedSum"} (how often a class borders anything). Ignored at
#'   patch scale.
#' @param connect [`integerish(1)`][integer]\cr which cells border each other:
#'   \code{4} (default) the orthogonal neighbours, \code{8} also the diagonal
#'   ones. At patch scale this decides only which patches touch; which cells
#'   form a patch was set by \code{\link{mdf_componentise}}.
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with the adjacencies added to the results of
#'   \code{layer}, as \code{likeAdj}, \code{pairedSum} or \code{adjacency} in
#'   its classes (see \code{\link{msk_categories}}), or as \code{adjacency}
#'   and \code{regions} in its patches (see \code{\link{msk_patches}}).
#' @details Each pair of bordering cells is counted from both sides, as in
#'   FRAGSTATS: two forest cells side by side add 2 to the like adjacencies of
#'   forest, a forest cell beside a meadow cell adds 1 to forest-meadow and 1
#'   to meadow-forest. The \code{"paired"} result is therefore a symmetric
#'   class by class matrix, its diagonal is \code{"like"} and its row sums are
#'   \code{"pairedSum"}.
#'
#'   At patch scale there are two patch by patch matrices, with the patch
#'   numbers as row and column names: \code{adjacency}, how many pairs of
#'   cells two patches share along their border, and \code{regions}, in how
#'   many separate places they touch. \code{regions[a, b]} counts the places
#'   from the side of patch \code{a}, so it need not equal \code{regions[b,
#'   a]}; a patch that touches another in two places forms a loop with it.
#' @examples
#' # how often each class borders itself
#' m <- msr_adjacency(landscape, type = "like")
#' msk_categories(m)$likeAdj
#'
#' # how often each class borders each other class
#' m <- msr_adjacency(landscape, type = "paired")
#' msk_categories(m)$adjacency
#'
#' # how often each class borders anything
#' m <- msr_adjacency(landscape, type = "pairedSum")
#' msk_categories(m)$pairedSum
#'
#' # to calculate patch-level metrics, number the patches first; scale =
#' # "patch" looks for the table mdf_componentise writes
#' f <- mdf_filter(landscape, cover == 47, add = "forest") |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
#'   msr_adjacency(scale = "patch", connect = 8, layer = "forest")
#' msk_patches(f, layer = "forest")$adjacency
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_adjacency <- function(obj = NULL, scale = "class", type = "like",
                          connect = 4, layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = scale, choices = c("class", "patch"))
  assertChoice(x = type, choices = c("like", "paired", "pairedSum"))
  assertChoice(x = connect, choices = c(4, 8))
  assertCharacter(x = layer, null.ok = TRUE)

  eightConn <- connect == 8

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims

  # patch scale: the patches numbered by mdf_componentise, then contact +
  # region matrices ----
  if(scale == "patch"){

    p <- .patches_of(obj, layer)
    labels <- p$ids
    labels[!labels %in% p$patch] <- NA

    pa <- patchAdjacencyCpp(labels = as.integer(labels), nrow = dims[2],
                           ncol = dims[1], eightconn = eightConn)

    ids <- pa$ids
    dimnames(pa$adjacency) <- list(as.character(ids), as.character(ids))
    dimnames(pa$regions)   <- list(as.character(ids), as.character(ids))

    # a patch the map border cuts may have contacts beyond it; a whole
    # patch's contacts are all on the map
    cut <- intersect(as.character(p$patch[p$clipped]), as.character(ids))
    pa$adjacency[cut, ] <- NA
    pa$regions[cut, ] <- NA

    obj@patches[[layer]]$adjacency <- pa$adjacency
    obj@patches[[layer]]$regions   <- pa$regions

    obj <- .update_mosaik(obj, step = step)
    return(obj)
  }

  uVals <- sort(unique(vals[!is.na(vals)]))

  # every pair of bordering cells counted from both sides
  values <- countCellAdjacenciesCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                                    doublecount = TRUE, eightconn = eightConn)
  rownames(values) <- uVals
  colnames(values) <- uVals

  # attach to the categories table, aligned with the classes it lists
  newGids <- as.integer(uVals)
  existing <- obj@categories[[layer]]
  if (is.null(existing)) existing <- list(gid = newGids)
  idx <- match(existing$gid, newGids)

  if(type == "like"){
    metric <- diag(values)
    existing$likeAdj <- ifelse(is.na(idx), NA, metric[idx])
  } else if(type == "paired"){
    # a class missing from the layer gets an NA row and column
    existing$adjacency <- values[idx, idx, drop = FALSE]
    dimnames(existing$adjacency) <- list(existing$gid, existing$gid)
  } else {
    metric <- rowSums(values)
    existing$pairedSum <- ifelse(is.na(idx), NA, metric[idx])
  }
  obj@categories[[layer]] <- existing

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
