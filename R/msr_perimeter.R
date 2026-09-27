#' Perimeter length
#'
#' Calculate the length of boundaries between objects in a mosaik and attach
#' the result to the attribute table.
#' @param obj [mosaik]\cr the mosaik to measure.
#' @param scale [character(1)][character]\cr scale at which to calculate;
#'   \code{"patch"}, \code{"class"} or \code{"landscape"}.
#' @param unit [character(1)][character]\cr \code{"cells"} (default) or
#'   \code{"map"} (in map units).
#' @param layer [character(1)][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with a \code{perimeter} column added to
#'   \code{@patches} (for \code{scale = "patch"}),
#'   \code{@categories} (for \code{scale = "class"}), or
#'   \code{@global} (for \code{scale = "landscape"}).
#' @examples
#' # landscape-level: total perimeter
#' m <- msr_perimeter(landscape, scale = "landscape")
#' m@global$perimeter
#'
#' # class-level: perimeter per class
#' m <- msr_perimeter(landscape, scale = "class")
#' m@categories$cover$perimeter
#'
#' # patch-level: perimeter per patch
#' m <- msr_perimeter(landscape, scale = "patch")
#' m@patches$perimeter
#'
#' # derive: edge density (perimeter / landscape area)
#' m <- msr_perimeter(landscape, scale = "class")
#' m <- msr_area(m, scale = "landscape")
#' m <- msr(m, equation = "perimeter.class / area.landscape",
#'          label = "edge_density")
#' m@categories$cover$edge_density
#'
#' # derive: shape index per patch
#' m <- msr_perimeter(landscape, scale = "patch")
#' m <- msr_area(m, scale = "patch")
#' m <- msr(m, equation = "perimeter.patch / sqrt(area.patch)",
#'          label = "shape_index")
#' m@patches$shape_index
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_perimeter <- function(obj, scale = "patch", unit = "cells", layer = NULL){

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

    obj@global$perimeter <- total

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

    # ensure cached patch layer exists (also sets @patches$class, $patch)
    obj <- .ensure_patch_layer(obj, layer)
    patch_ids <- msk_pull(obj, "_patches")
    patch_ids[is.na(patch_ids)] <- 0

    # count edges per patch (globally unique IDs)
    temp_edges <- countCellEdgesCpp(vals = as.numeric(patch_ids),
                                    nrow = dims[1], ncol = dims[2])
    temp_edges <- temp_edges[temp_edges$value != 0,]

    # align with the roster order from .ensure_patch_layer
    idx <- match(obj@patches$patch, temp_edges$value)
    if(unit == "map"){
      obj@patches$perimeter <- temp_edges$edgesX[idx] * theRes[1] +
                               temp_edges$edgesY[idx] * theRes[2]
    } else {
      obj@patches$perimeter <- temp_edges$edgesX[idx] + temp_edges$edgesY[idx]
    }

  }

  # provenance
  prov <- msk_prov("msr_perimeter", list(scale = scale, unit = unit, layer = layer))
  obj@provenance <- c(obj@provenance, list(prov))

  return(obj)
}
