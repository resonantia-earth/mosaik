# TODO: hm, something might be wrong here, the "percolation" doesn't really look
# like clusters. Just a mixed up, slightly clustered "noise".

#' Make a clusters surface
#'
#' Make categorical landscapes via stochastic spatial growth processes.
#' \code{"percolation"} creates a binary random map at a given probability.
#' \code{"randomCluster"} extends percolation with connected-component
#' labelling and class assignment.
#'
#' @param obj [mosaik]\cr the mosaik whose grid the clusters are generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [character(1)][character]\cr the cluster type:
#'   \code{"percolation"} (default) or \code{"randomCluster"}.
#' @param p [numeric(1)][numeric]\cr percolation probability (proportion of
#'   cells set to foreground). Default 0.5.
#' @param n [integerish(1)][integer]\cr number of classes for
#'   \code{"randomCluster"}. Default 3. Ignored for \code{"percolation"}.
#' @param name [character(1)][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param role [character(1)][character]\cr optional semantic role tag (e.g.
#'   \code{"surface"}, \code{"precipitation"}). Used by downstream functions to
#'   auto-detect layer purpose.
#' @param seed [integerish(1)][integer]\cr random seed for reproducibility.
#' @return A \code{\link{mosaik}} with integer cluster labels.
#' @references
#'   Stauffer D, Aharony A. Introduction to Percolation Theory. 2nd ed.
#'   Taylor & Francis; 1994.
#'
#'   Saura S, Martinez-Millan J. Landscape patterns simulation with a modified
#'   random clusters method. Landscape Ecology. 2000;15:661-678.
#' @family generator functions
#' @examples
#' m <- mosaik(extent = c(0, 100, 0, 100), res = 1)
#'
#' # binary percolation at the critical threshold
#' syn_cluster(m, type = "percolation", p = 0.5, seed = 1)
#'
#' # random cluster landscape with 4 land-cover classes
#' syn_cluster(m, type = "randomCluster", p = 0.6, n = 4, seed = 1)
#' @importFrom checkmate assertClass assertIntegerish assertChoice assertNumber
#'   assertCharacter
#' @export

syn_cluster <- function(obj, type = "percolation", p = 0.5, n = 3L,
                         name = "values", role = NULL, seed = NULL) {

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type, choices = c("percolation", "randomCluster"))
  assertNumber(x = p, lower = 0, upper = 1)
  assertIntegerish(x = n, len = 1, lower = 2)
  assertCharacter(x = name, len = 1)
  assertCharacter(x = role, len = 1, null.ok = TRUE)
  assertIntegerish(x = seed, len = 1, null.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  ncols <- obj@dims[1]
  nrows <- obj@dims[2]
  n_cells <- ncols * nrows

  vals <- switch(type,

    percolation = {
      as.numeric(percolationCpp(ncol = ncols, nrow = nrows, p = p))
    },

    randomCluster = {
      binary <- percolationCpp(ncol = ncols, nrow = nrows, p = p)

      comps <- componentsCpp(vals = as.numeric(binary),
                             nrow = nrows, ncol = ncols)

      ucomps <- sort(unique(comps[comps > 0]))
      class_map <- sample.int(as.integer(n), length(ucomps), replace = TRUE)
      names(class_map) <- as.character(ucomps)

      out <- rep(0L, n_cells)
      for (i in seq_along(ucomps)) {
        out[comps == ucomps[i]] <- class_map[i]
      }

      if (any(out == 0L)) {
        filled <- out
        while (any(filled == 0L)) {
          zeros <- which(filled == 0L)
          for (idx in zeros) {
            row_i <- ((idx - 1L) %/% ncols)
            col_i <- ((idx - 1L) %% ncols)
            neighbours <- integer(0)
            if (col_i > 0)         neighbours <- c(neighbours, filled[idx - 1L])
            if (col_i < ncols - 1) neighbours <- c(neighbours, filled[idx + 1L])
            if (row_i > 0)         neighbours <- c(neighbours, filled[idx - ncols])
            if (row_i < nrows - 1) neighbours <- c(neighbours, filled[idx + ncols])
            nn <- neighbours[neighbours > 0L]
            if (length(nn) > 0) {
              filled[idx] <- nn[sample.int(length(nn), 1)]
            }
          }
        }
        out <- filled
      }

      as.numeric(out)
    }
  )

  prov <- msk_prov("syn_cluster",
                     list(type = type, p = p, n = n, name = name, seed = seed))

  msk_set(obj, name, vals, prov, role = role)
}
