#' Area of objects
#'
#' Calculate the area of objects in a mosaik and attach the result to the
#' attribute table.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param scale [`character(1)`][character]\cr scale at which to calculate;
#'   \code{"patch"}, \code{"class"} or \code{"landscape"}.
#' @param unit [`character(1)`][character]\cr \code{"cells"} (default, number of
#'   cells) or \code{"map"} (in map units).
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with \code{area} added to the results of
#'   \code{layer}, in its patches, classes or landscape values (see
#'   \code{\link{msk_patches}}, \code{\link{msk_categories}},
#'   \code{\link{msk_global}}).
#' @details At \code{scale = "patch"}, the patches are those numbered by
#'   \code{\link{mdf_componentise}} on \code{layer}, which must run first.
#' @examples
#' # the area of the whole layer
#' m <- msr_area(landscape, scale = "landscape")
#' msk_global(m)$area
#'
#' # the area of each class
#' m <- msr_area(landscape, scale = "class")
#' msk_categories(m)$area
#'
#' # to calculate patch-level metrics, number the patches first; scale =
#' # "patch" looks for the table mdf_componentise writes
#' m <- mdf_componentise(landscape, connectivity = 8L, layer = "cover",
#'                       add = "patch")
#' m <- msr_area(m, scale = "patch")
#' msk_patches(m)$area
#' msk_patches(m)$class   # which class each patch belongs to
#' msk_patches(m)$patch   # the patch number, as in the layer 'patch'
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_area <- function(obj = NULL, scale = "patch", unit = "cells", layer = NULL){

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

    total <- sum(!is.na(vals))
    if(unit == "map"){
      total <- total * theRes[1] * theRes[2]
    }

    # attach to the layer's landscape-level values
    obj@global[[layer]]$area <- total

  } else if(scale == "class"){

    values <- countCellValuesCpp(vals = vals, nrow = dims[2], ncol = dims[1])
    classes <- values$value

    if(unit == "map"){
      area <- values$cells * theRes[1] * theRes[2]
    } else{
      area <- values$cells
    }

    # attach to categories table
    newGids <- as.integer(classes)
    existing <- obj@categories[[layer]]
    if (!is.null(existing)) {
      idx <- match(existing$gid, newGids)
      existing$area <- ifelse(is.na(idx), NA, area[idx])
      obj@categories[[layer]] <- existing
    } else {
      obj@categories[[layer]] <- list(gid = newGids, area = area)
    }

  } else {

    # the patches numbered by mdf_componentise
    p <- .patches_of(obj, layer)
    patch_ids <- p$ids
    patch_ids[!patch_ids %in% p$patch] <- NA

    # count cells per patch
    temp_counts <- countCellValuesCpp(vals = as.numeric(patch_ids),
                                      nrow = dims[2], ncol = dims[1])
    temp_counts <- temp_counts[!is.na(temp_counts$value),]

    # align with the order of the patch record
    idx <- match(p$patch, temp_counts$value)
    if(unit == "map"){
      area <- temp_counts$cells[idx] * theRes[1] * theRes[2]
    } else {
      area <- temp_counts$cells[idx]
    }
    area[p$clipped] <- NA
    obj@patches[[layer]]$area <- area

  }

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
