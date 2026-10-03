#' Make a clusters surface
#'
#' Make categorical landscapes from percolation. \code{"percolation"} sets
#' every cell to 1 with probability \code{p}, independently of all others;
#' the clusters are the connected areas of 1s, which grow with \code{p} and
#' first span the whole map at about \code{p = 0.593} (with four neighbours).
#' The map itself therefore looks like noise, and its clusters show once they
#' are numbered with \code{\link{mdf_componentise}}. \code{"randomCluster"}
#' (the modified random clusters method) gives each cluster one of \code{n}
#' classes at random and fills the cells between the clusters from their
#' neighbours, which makes a categorical map of irregular patches.
#'
#' @param obj [`mosaik`]\cr the mosaik whose grid the clusters are generated on.
#'   Create an empty grid with \code{\link{mosaik}(extent, res)}.
#' @param type [`character(1)`][character]\cr the cluster type:
#'   \code{"percolation"} (default) or \code{"randomCluster"}.
#' @param p [`numeric(1)`][numeric]\cr percolation probability (proportion of
#'   cells set to foreground). Default 0.5.
#' @param n [`integerish(1)`][integer]\cr number of classes for
#'   \code{"randomCluster"}. Default 3. Ignored for \code{"percolation"}.
#' @param name [`character(1)`][character]\cr the layer name. Default
#'   \code{"values"}.
#' @param seed [`integerish(1)`][integer]\cr random seed for reproducibility.
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
#' # percolation near the threshold, and its clusters numbered (drawn in
#' # shuffled colours, since neighbouring numbers would get similar ones)
#' cl <- m |>
#'   drw_cluster(type = "percolation", p = 0.59, seed = 1,
#'               name = "percolation") |>
#'   mdf_componentise(layer = "percolation", add = "clusters")
#' set.seed(1)
#' colours <- sample(hcl.colors(max(msk_pull(cl, "clusters"), na.rm = TRUE),
#'                              "Dark 3"))
#' msk_vis(cl, .layer("percolation"),
#'         .layer("clusters", colours = colours, legend = FALSE))
#'
#' # random cluster maps with 4 classes, from sparser and denser percolation
#' rc <- m |>
#'   drw_cluster(type = "randomCluster", p = 0.4, n = 4, seed = 1,
#'               name = "p_0.4") |>
#'   drw_cluster(type = "randomCluster", p = 0.55, n = 4, seed = 1,
#'               name = "p_0.55")
#' msk_vis(rc, .layer("p_0.4"), .layer("p_0.55"))
#' @importFrom checkmate assertClass assertIntegerish assertChoice assertNumber
#'   assertCharacter
#' @export

drw_cluster <- function(obj, type = "percolation", p = 0.5, n = 3L,
                         name = "values", seed = NULL) {

  step <- .step()

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type, choices = c("percolation", "randomCluster"))
  assertNumber(x = p, lower = 0, upper = 1)
  assertIntegerish(x = n, len = 1, lower = 2)
  assertCharacter(x = name, len = 1)
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

      # only the cells that percolated form clusters; the others are filled
      # from their neighbours below. Numbering the 0 cells as well would give
      # the gaps between the clusters classes of their own
      fg <- as.numeric(binary)
      fg[fg == 0] <- NA
      comps <- componentsCpp(vals = fg, nrow = nrows, ncol = ncols)
      comps[is.na(comps)] <- 0L

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


  .update_mosaik(obj, values = vals, step = step)
}
