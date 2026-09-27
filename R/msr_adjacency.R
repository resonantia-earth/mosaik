#' Adjacency of cells
#'
#' Calculate the cell-adjacency matrix for a mosaik and attach the result to
#' the attribute table.
#' @param obj [mosaik]\cr the mosaik to measure.
#' @param scale [character(1)][character]\cr \code{"class"} (default) measures
#'   adjacency between class values; \code{"patch"} labels every patch across
#'   the whole grid as a distinct node (global componentisation) and measures
#'   adjacency between patches. \code{type} is ignored at patch scale.
#' @param type [character(1)][character]\cr (class scale only) which adjacencies
#'   to calculate; \code{"like"} (diagonal of adjacency matrix), \code{"paired"}
#'   (full matrix) or \code{"pairedSum"} (row sums).
#' @param count [character(1)][character]\cr \code{"single"} counts only right
#'   and bottom neighbours; \code{"double"} also counts left and top.
#' @param connect [integerish(1)][integer]\cr neighbourhood rule: \code{4}
#'   (rook, orthogonal neighbours only) or \code{8} (queen, also diagonal
#'   neighbours). Defaults to \code{4}. At patch scale this governs both the
#'   componentisation and what counts as contact between patches.
#' @param layer [character(1)][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with adjacency values attached. At \code{scale =
#'   "class"}: \code{@categories} (for \code{type = "like"}/\code{"pairedSum"})
#'   or \code{@global} (for \code{type = "paired"}, as a matrix). At \code{scale
#'   = "patch"}: two patch \eqn{\times} patch matrices in \code{@patches} —
#'   \code{adjacency} (count of adjacent cell pairs between patches, symmetric)
#'   and \code{regions} (number of spatially distinct contact places, read
#'   \code{regions[fragment, core]}; asymmetric). Row/column names are the
#'   global patch IDs.
#' @examples
#' # like-adjacency: number of same-class cell pairs per class
#' m <- msr_adjacency(landscape, type = "like")
#' m@categories$cover$likeAdj
#'
#' # full adjacency matrix (landscape-level)
#' m <- msr_adjacency(landscape, type = "paired")
#' m@global$adjacency
#'
#' # row sums of adjacency matrix (per class)
#' m <- msr_adjacency(landscape, type = "pairedSum")
#' m@categories$cover$pairedSum
#'
#' # derive: percentage of like adjacencies per class
#' m <- msr_adjacency(landscape, type = "like")
#' m <- msr_adjacency(m, type = "pairedSum")
#' m <- msr(m, equation = "likeAdj.class / pairedSum.class * 100",
#'          label = "pladj")
#' m@categories$cover$pladj
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_adjacency <- function(obj, scale = "class", type = "like", count = "double",
                          connect = 4, layer = NULL){

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

  # patch scale: global component labelling, then contact + region matrices ----
  if(scale == "patch"){

    labels <- componentsCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                            connectivity = connect)

    pa <- patchAdjacencyCpp(labels = as.integer(labels), nrow = dims[2],
                           ncol = dims[1], eightconn = eightConn)

    ids <- pa$ids
    dimnames(pa$adjacency) <- list(as.character(ids), as.character(ids))
    dimnames(pa$regions)   <- list(as.character(ids), as.character(ids))

    obj@patches$adjacency <- pa$adjacency
    obj@patches$regions   <- pa$regions

    prov <- msk_prov("msr_adjacency", list(scale = scale, connect = connect,
                       layer = layer))
    obj@provenance <- c(obj@provenance, list(prov))
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
    obj@global$adjacency <- values

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
  prov <- msk_prov("msr_adjacency", list(type = type, count = count,
                     connect = connect, layer = layer))
  obj@provenance <- c(obj@provenance, list(prov))

  return(obj)
}
