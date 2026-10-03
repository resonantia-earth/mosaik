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
#' @return The input mosaik with an \code{area} value added to the results of
#'   \code{layer}: its patches (for \code{scale = "patch"}, which must have
#'   been numbered with \code{\link{mdf_componentise}} first, see
#'   \code{\link{msk_patches}}), its classes (\code{scale = "class"}, see
#'   \code{\link{msk_categories}}) or the layer as a whole
#'   (\code{scale = "landscape"}, see \code{\link{msk_global}}).
#' @examples
#' # landscape-level: total area
#' m <- msr_area(landscape, scale = "landscape")
#' msk_global(m)$area
#'
#' # class-level: area per class
#' m <- msr_area(landscape, scale = "class")
#' msk_categories(m)$area
#'
#' # patch-level: area per patch, of the patches numbered first
#' m <- mdf_componentise(landscape, connectivity = 8L, layer = "cover",
#'                       add = "patch")
#' m <- msr_area(m, scale = "patch")
#' msk_patches(m)$area
#' msk_patches(m)$class   # which class each patch belongs to
#' msk_patches(m)$patch   # the patch number, as in the layer 'patch'
#'
#' # derive: proportion of landscape per class (PLAND)
#' m <- msr_area(m, scale = "class")
#' m <- msr_area(m, scale = "landscape")
#' m <- msr(m, equation = "area.class / area.landscape * 100",
#'          label = "pland")
#' msk_categories(m)$pland
#'
#' # derive: mean patch size per class
#' m <- msr_number(m, scale = "patch")
#' m <- msr(m, equation = "area.class / number.class",
#'          label = "mean_patch_size")
#' msk_categories(m)$mean_patch_size
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

    values <- countCellValuesCpp(vals = vals, nrow = dims[1], ncol = dims[2])
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
                                      nrow = dims[1], ncol = dims[2])
    temp_counts <- temp_counts[!is.na(temp_counts$value),]

    # align with the order of the patch record
    idx <- match(p$patch, temp_counts$value)
    if(unit == "map"){
      obj@patches[[layer]]$area <- temp_counts$cells[idx] * theRes[1] * theRes[2]
    } else {
      obj@patches[[layer]]$area <- temp_counts$cells[idx]
    }

  }

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
