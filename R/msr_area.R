#' Area of objects
#'
#' Calculate the area of objects in a mosaik and attach the result to the
#' attribute table.
#' @param obj [mosaik]\cr the mosaik to measure.
#' @param scale [character(1)][character]\cr scale at which to calculate;
#'   \code{"patch"}, \code{"class"} or \code{"landscape"}.
#' @param unit [character(1)][character]\cr \code{"cells"} (default, number of
#'   cells) or \code{"map"} (in map units).
#' @param layer [character(1)][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with an \code{area} column added to
#'   \code{@patches} (for \code{scale = "patch"}),
#'   \code{@categories} (for \code{scale = "class"}), or
#'   \code{@global} (for \code{scale = "landscape"}).
#' @examples
#' # landscape-level: total area
#' m <- msr_area(landscape, scale = "landscape")
#' m@global$area
#'
#' # class-level: area per class
#' m <- msr_area(landscape, scale = "class")
#' m@categories$cover$area
#'
#' # patch-level: area per patch
#' m <- msr_area(landscape, scale = "patch")
#' m@patches$area
#' m@patches$class   # which class each patch belongs to
#' m@patches$patch    # patch ID within each class
#'
#' # derive: proportion of landscape per class (PLAND)
#' m <- msr_area(landscape, scale = "class")
#' m <- msr_area(m, scale = "landscape")
#' m <- msr(m, equation = "area.class / area.landscape * 100",
#'          label = "pland")
#' m@categories$cover$pland
#'
#' # derive: mean patch size per class
#' m <- msr_number(m, scale = "patch")
#' m <- msr(m, equation = "area.class / number.class",
#'          label = "mean_patch_size")
#' m@categories$cover$mean_patch_size
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_area <- function(obj, scale = "patch", unit = "cells", layer = NULL){

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

    # attach to global table
    obj@global$area <- total

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

    # ensure cached patch layer exists (also sets @patches$class, $patch)
    obj <- .ensure_patch_layer(obj, layer)
    patch_ids <- msk_pull(obj, "_patches")

    # count cells per patch
    temp_counts <- countCellValuesCpp(vals = as.numeric(patch_ids),
                                      nrow = dims[1], ncol = dims[2])
    temp_counts <- temp_counts[!is.na(temp_counts$value),]

    # align with the roster order from .ensure_patch_layer
    idx <- match(obj@patches$patch, temp_counts$value)
    if(unit == "map"){
      obj@patches$area <- temp_counts$cells[idx] * theRes[1] * theRes[2]
    } else {
      obj@patches$area <- temp_counts$cells[idx]
    }

  }

  # provenance
  prov <- msk_prov("msr_area", list(scale = scale, unit = unit, layer = layer))
  obj@provenance <- c(obj@provenance, list(prov))

  return(obj)
}
