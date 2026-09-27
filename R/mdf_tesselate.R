#' Create a Voronoi tessellation from patches
#'
#' Assign every cell to the nearest patch centroid, creating a Voronoi-like
#' partitioning of the grid. The input must contain distinct patches (run
#' \code{\link{mdf_componentise}} first if needed).
#' @param obj [mosaik]\cr the mosaik to tessellate. Non-NA foreground cells
#'   are treated as patches.
#' @param layer [character(1)][character]\cr the layer to use. Defaults to the
#'   first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions where every cell carries the ID of
#'   the nearest patch centroid.
#' @examples
#' forest <- mdf_binarise(landscape, match = 47, layer = "cover")
#' patches <- mdf_componentise(forest)
#' mdf_tesselate(patches)
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertCharacter
#' @export

mdf_tesselate <- function(obj = NULL,
                          layer = NULL,
                          add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  ext  <- obj@extent
  res  <- msk_res(obj)

  # find patch centroids ----
  uvals <- sort(unique(vals[!is.na(vals)]))
  if(length(uvals) < 2){
    stop("Need at least 2 distinct patch values to tessellate.")
  }

  xs <- seq(ext[1] + res[1] / 2, by = res[1], length.out = dims[1])
  ys <- seq(ext[4] - res[2] / 2, by = -res[2], length.out = dims[2])
  coords_x <- rep(xs, times = dims[2])
  coords_y <- rep(ys, each = dims[1])

  # compute centroid for each patch value
  seed_x <- numeric(length(uvals))
  seed_y <- numeric(length(uvals))
  seed_val <- as.integer(uvals)
  for(i in seq_along(uvals)){
    mask <- which(vals == uvals[i])
    seed_x[i] <- mean(coords_x[mask])
    seed_y[i] <- mean(coords_y[mask])
  }

  # assign each cell to nearest centroid (C++) ----
  temp <- tesselateCpp(xmin = ext[1], ymax = ext[4],
                       res_x = res[1], res_y = res[2],
                       ncol = dims[1], nrow = dims[2],
                       seed_x = seed_x, seed_y = seed_y,
                       seed_val = seed_val)

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_tesselate", list(n_patches = length(uvals),
                                            layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
