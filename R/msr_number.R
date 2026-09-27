#' Number of objects
#'
#' Count the number of objects at the given scale in a mosaik and attach the
#' result to the attribute table.
#' @param obj [mosaik]\cr the mosaik to measure.
#' @param scale [character(1)][character]\cr the level to report at;
#'   \code{"class"} (number of distinct classes, attached to
#'   \code{@global}) or \code{"patch"} (number of patches per class,
#'   attached to \code{@categories}).
#' @param layer [character(1)][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with results attached to
#'   \code{@global} (for \code{scale = "class"}: the count of distinct
#'   classes) or \code{@categories} (for \code{scale = "patch"}: patch
#'   count per class).
#' @examples
#' # count distinct classes (landscape-level)
#' m <- msr_number(landscape, scale = "class")
#' m@global$number
#'
#' # count patches per class (class-level)
#' m <- msr_number(landscape, scale = "patch")
#' m@categories$cover$number
#'
#' # derive: patch density (patches per unit area)
#' m <- msr_number(landscape, scale = "patch")
#' m <- msr_area(m, scale = "landscape")
#' m <- msr(m, equation = "number.class / area.landscape",
#'          label = "patch_density")
#' m@categories$cover$patch_density
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_number <- function(obj, scale = "class", layer = NULL){

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = scale, choices = c("class", "patch"))
  assertCharacter(x = layer, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims

  uVals <- sort(unique(vals[!is.na(vals)]))

  if(scale == "class"){

    # number of distinct classes (landscape-level result)
    obj@global$number <- length(uVals)

  } else {

    # number of patches per class (class-level result)
    # reuse cached patch layer
    obj <- .ensure_patch_layer(obj, layer)
    patches_per_class <- as.integer(table(obj@patches$class))

    # attach patch count per class to categories
    newGids <- as.integer(uVals)
    existing <- obj@categories[[layer]]
    if (!is.null(existing)) {
      idx <- match(existing$gid, newGids)
      existing$number <- ifelse(is.na(idx), NA_integer_, patches_per_class[idx])
      obj@categories[[layer]] <- existing
    } else {
      obj@categories[[layer]] <- list(gid = newGids, number = patches_per_class)
    }

  }

  # provenance
  prov <- msk_prov("msr_number", list(scale = scale, layer = layer))
  obj@provenance <- c(obj@provenance, list(prov))

  return(obj)
}
