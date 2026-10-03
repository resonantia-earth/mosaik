#' Number of objects
#'
#' Count the number of objects at the given scale in a mosaik and attach the
#' result to the attribute table.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param scale [`character(1)`][character]\cr the level to report at;
#'   \code{"class"} (number of distinct classes, a landscape-level value) or
#'   \code{"patch"} (number of patches per class, a class-level value; the
#'   patches must have been numbered with \code{\link{mdf_componentise}}).
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with a \code{number} value added to the results of
#'   \code{layer}: at landscape level for \code{scale = "class"} (see
#'   \code{\link{msk_global}}), at class level for \code{scale = "patch"} (see
#'   \code{\link{msk_categories}}).
#' @examples
#' # count distinct classes (landscape-level)
#' m <- msr_number(landscape, scale = "class")
#' msk_global(m)$number
#'
#' # count patches per class (class-level), of the patches numbered first
#' p <- mdf_componentise(landscape, connectivity = 8L, layer = "cover",
#'                       add = "patch")
#' m <- msr_number(p, scale = "patch")
#' msk_categories(m)$number
#'
#' # derive: patch density (patches per unit area)
#' m <- msr_number(p, scale = "patch")
#' m <- msr_area(m, scale = "landscape")
#' m <- msr(m, equation = "number.class / area.landscape",
#'          label = "patch_density")
#' msk_categories(m)$patch_density
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_number <- function(obj = NULL, scale = "class", layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

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
    obj@global[[layer]]$number <- length(uVals)

  } else {

    # number of patches per class (class-level result), from the patches
    # numbered by mdf_componentise; a class without patches (0) counts none
    p <- .patches_of(obj, layer)
    patches_per_class <- vapply(uVals, function(v) sum(p$class == v),
                                integer(1))

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
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
