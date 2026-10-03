#' Adjacency of cells
#'
#' Calculate the cell-adjacency matrix for a mosaik and attach the result to
#' the attribute table.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param scale [`character(1)`][character]\cr \code{"class"} (default) measures
#'   adjacency between class values; \code{"patch"} measures adjacency between
#'   the patches of \code{layer}, which must have been numbered with
#'   \code{\link{mdf_componentise}} first. \code{type} is ignored at patch
#'   scale.
#' @param type [`character(1)`][character]\cr (class scale only) which adjacencies
#'   to calculate; \code{"like"} (diagonal of adjacency matrix), \code{"paired"}
#'   (full matrix) or \code{"pairedSum"} (row sums).
#' @param count [`character(1)`][character]\cr \code{"single"} counts only right
#'   and bottom neighbours; \code{"double"} also counts left and top.
#' @param connect [`integerish(1)`][integer]\cr neighbourhood rule: \code{4}
#'   (rook, orthogonal neighbours only) or \code{8} (queen, also diagonal
#'   neighbours). Defaults to \code{4}. At patch scale this governs what counts
#'   as contact between patches; which cells form a patch was set by
#'   \code{\link{mdf_componentise}}.
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with adjacency values added to the results of
#'   \code{layer}. At \code{scale = "class"}: class level (for \code{type =
#'   "like"}/\code{"pairedSum"}, see \code{\link{msk_categories}}) or landscape
#'   level (for \code{type = "paired"}, as a matrix, see
#'   \code{\link{msk_global}}). At \code{scale = "patch"}: two patch
#'   \eqn{\times} patch matrices at patch level (see \code{\link{msk_patches}}) —
#'   \code{adjacency} (count of adjacent cell pairs between patches, symmetric)
#'   and \code{regions} (number of spatially distinct contact places, read
#'   \code{regions[fragment, core]}; asymmetric). Row/column names are the
#'   patch numbers.
#' @examples
#' # like-adjacency: number of same-class cell pairs per class
#' m <- msr_adjacency(landscape, type = "like")
#' msk_categories(m)$likeAdj
#'
#' # full adjacency matrix (landscape-level)
#' m <- msr_adjacency(landscape, type = "paired")
#' msk_global(m)$adjacency
#'
#' # row sums of adjacency matrix (per class)
#' m <- msr_adjacency(landscape, type = "pairedSum")
#' msk_categories(m)$pairedSum
#'
#' # derive: percentage of like adjacencies per class
#' m <- msr_adjacency(landscape, type = "like")
#' m <- msr_adjacency(m, type = "pairedSum")
#' m <- msr(m, equation = "likeAdj.class / pairedSum.class * 100",
#'          label = "pladj")
#' msk_categories(m)$pladj
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_adjacency <- function(obj = NULL, scale = "class", type = "like", count = "double",
                          connect = 4, layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = scale, choices = c("class", "patch"))
  assertChoice(x = type, choices = c("like", "paired", "pairedSum"))
  assertChoice(x = count, choices = c("single", "double"))
  assertChoice(x = connect, choices = c(4, 8))
  assertCharacter(x = layer, null.ok = TRUE)

  countDouble <- count == "double"
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

    obj@patches[[layer]]$adjacency <- pa$adjacency
    obj@patches[[layer]]$regions   <- pa$regions

    obj <- .update_mosaik(obj, step = step)
    return(obj)
  }

  uVals <- sort(unique(vals[!is.na(vals)]))

  values <- countCellAdjacenciesCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                                    doublecount = countDouble, eightconn = eightConn)
  rownames(values) <- uVals
  colnames(values) <- uVals

  if(type == "like"){

    # attach to categories table
    newGids <- as.integer(uVals)
    metric <- diag(values)
    existing <- obj@categories[[layer]]
    if (!is.null(existing)) {
      idx <- match(existing$gid, newGids)
      existing$likeAdj <- ifelse(is.na(idx), NA, metric[idx])
      obj@categories[[layer]] <- existing
    } else {
      obj@categories[[layer]] <- list(gid = newGids, likeAdj = metric)
    }

  } else if(type == "paired"){

    # attach full adjacency matrix to global table
    obj@global[[layer]]$adjacency <- values

  } else {

    # attach to categories table
    newGids <- as.integer(uVals)
    metric <- rowSums(values)
    existing <- obj@categories[[layer]]
    if (!is.null(existing)) {
      idx <- match(existing$gid, newGids)
      existing$pairedSum <- ifelse(is.na(idx), NA, metric[idx])
      obj@categories[[layer]] <- existing
    } else {
      obj@categories[[layer]] <- list(gid = newGids, pairedSum = metric)
    }

  }

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
