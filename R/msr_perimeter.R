#' Perimeter length
#'
#' Calculate the length of boundaries between objects in a mosaik and attach
#' the result to the attribute table.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param scale [`character(1)`][character]\cr scale at which to calculate;
#'   \code{"patch"}, \code{"class"} or \code{"landscape"}.
#' @param unit [`character(1)`][character]\cr \code{"cells"} (default) or
#'   \code{"map"} (in map units).
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with a \code{perimeter} value added to the results
#'   of \code{layer}: its patches (for \code{scale = "patch"}, which must have
#'   been numbered with \code{\link{mdf_componentise}} first, see
#'   \code{\link{msk_patches}}), its classes (\code{scale = "class"}, see
#'   \code{\link{msk_categories}}) or the layer as a whole
#'   (\code{scale = "landscape"}, see \code{\link{msk_global}}).
#' @examples
#' # landscape-level: total perimeter
#' m <- msr_perimeter(landscape, scale = "landscape")
#' msk_global(m)$perimeter
#'
#' # class-level: perimeter per class
#' m <- msr_perimeter(landscape, scale = "class")
#' msk_categories(m)$perimeter
#'
#' # patch-level: perimeter per patch, of the patches numbered first
#' p <- mdf_componentise(landscape, connectivity = 8L, layer = "cover",
#'                       add = "patch")
#' m <- msr_perimeter(p, scale = "patch")
#' msk_patches(m)$perimeter
#'
#' # derive: edge density (perimeter / landscape area)
#' m <- msr_perimeter(landscape, scale = "class")
#' m <- msr_area(m, scale = "landscape")
#' m <- msr(m, equation = "perimeter.class / area.landscape",
#'          label = "edge_density")
#' msk_categories(m)$edge_density
#'
#' # derive: shape index per patch
#' m <- msr_perimeter(p, scale = "patch")
#' m <- msr_area(m, scale = "patch")
#' m <- msr(m, equation = "perimeter.patch / sqrt(area.patch)",
#'          label = "shape_index")
#' msk_patches(m)$shape_index
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_perimeter <- function(obj = NULL, scale = "patch", unit = "cells", layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = scale, choices = c("patch", "class", "landscape"))
  assertChoice(x = unit, choices = c("cells", "map"))
  assertCharacter(x = layer, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  theRes <- msk_res(obj)

  if(scale == "landscape"){

    values <- countCellEdgesCpp(vals = vals, nrow = dims[1], ncol = dims[2])
    if(unit == "map"){
      total <- sum(values$edgesX) * theRes[1] + sum(values$edgesY) * theRes[2]
    } else {
      total <- sum(values$edgesX) + sum(values$edgesY)
    }

    obj@global[[layer]]$perimeter <- total

  } else if(scale == "class"){

    values <- countCellEdgesCpp(vals = vals, nrow = dims[1], ncol = dims[2])
    classes <- values$value

    if(unit == "map"){
      edges <- values$edgesX * theRes[1] + values$edgesY * theRes[2]
    } else {
      edges <- values$edgesX + values$edgesY
    }

    # attach to categories table
    newGids <- as.integer(classes)
    existing <- obj@categories[[layer]]
    if (!is.null(existing)) {
      idx <- match(existing$gid, newGids)
      existing$perimeter <- ifelse(is.na(idx), NA, edges[idx])
      obj@categories[[layer]] <- existing
    } else {
      obj@categories[[layer]] <- list(gid = newGids, perimeter = edges)
    }

  } else {

    # the patches numbered by mdf_componentise; cells outside every patch
    # (including a 'background' value) become 0, which is not a patch number
    p <- .patches_of(obj, layer)
    patch_ids <- p$ids
    patch_ids[is.na(patch_ids) | !patch_ids %in% p$patch] <- 0

    # count edges per patch
    temp_edges <- countCellEdgesCpp(vals = as.numeric(patch_ids),
                                    nrow = dims[1], ncol = dims[2])
    temp_edges <- temp_edges[temp_edges$value != 0,]

    # align with the order of the patch record
    idx <- match(p$patch, temp_edges$value)
    if(unit == "map"){
      obj@patches[[layer]]$perimeter <- temp_edges$edgesX[idx] * theRes[1] +
                                        temp_edges$edgesY[idx] * theRes[2]
    } else {
      obj@patches[[layer]]$perimeter <- temp_edges$edgesX[idx] + temp_edges$edgesY[idx]
    }

  }

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
