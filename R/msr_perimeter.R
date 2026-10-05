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
#' @return The input mosaik with \code{perimeter} added to the results of
#'   \code{layer}, in its patches, classes or landscape values (see
#'   \code{\link{msk_patches}}, \code{\link{msk_categories}},
#'   \code{\link{msk_global}}).
#' @details At \code{scale = "patch"}, the patches are those numbered by
#'   \code{\link{mdf_componentise}} on \code{layer}, which must run first.
#'   A patch that touches the map border has perimeter \code{NA}: the map
#'   ends there, the patch may not, so its outline is unknown.
#' @examples
#' # the perimeter of the whole layer
#' m <- msr_perimeter(landscape, scale = "landscape")
#' msk_global(m)$perimeter
#'
#' # the perimeter of each class
#' m <- msr_perimeter(landscape, scale = "class")
#' msk_categories(m)$perimeter
#'
#' # to calculate patch-level metrics, number the patches first; scale =
#' # "patch" looks for the table mdf_componentise writes
#' p <- mdf_componentise(landscape, connectivity = 8L, layer = "cover",
#'                       add = "patch")
#' m <- msr_perimeter(p, scale = "patch")
#' msk_patches(m)$perimeter
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

    values <- countCellEdgesCpp(vals = vals, nrow = dims[2], ncol = dims[1])
    if(unit == "map"){
      total <- sum(values$edgesX) * theRes[2] + sum(values$edgesY) * theRes[1]
    } else {
      total <- sum(values$edgesX) + sum(values$edgesY)
    }

    obj@global[[layer]]$perimeter <- total

  } else if(scale == "class"){

    values <- countCellEdgesCpp(vals = vals, nrow = dims[2], ncol = dims[1])
    classes <- values$value

    if(unit == "map"){
      edges <- values$edgesX * theRes[2] + values$edgesY * theRes[1]
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
                                    nrow = dims[2], ncol = dims[1])
    temp_edges <- temp_edges[temp_edges$value != 0,]

    # align with the order of the patch record
    idx <- match(p$patch, temp_edges$value)
    if(unit == "map"){
      perim <- temp_edges$edgesX[idx] * theRes[2] + temp_edges$edgesY[idx] * theRes[1]
    } else {
      perim <- temp_edges$edgesX[idx] + temp_edges$edgesY[idx]
    }

    perim[p$clipped] <- NA
    obj@patches[[layer]]$perimeter <- perim

  }

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
