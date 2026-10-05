#' Number of objects
#'
#' Count the classes of a layer, or the patches of each class.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param scale [`character(1)`][character]\cr where the count is stored:
#'   \code{"landscape"} (default) counts the classes of the layer,
#'   \code{"class"} counts the patches of each class, which must have been
#'   numbered with \code{\link{mdf_componentise}} first.
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with \code{number} added to the results of
#'   \code{layer}, in its landscape values (see \code{\link{msk_global}}) or
#'   in its classes (see \code{\link{msk_categories}}).
#' @details A class that forms no patch, such as 0 on a binary layer, counts
#'   0 patches.
#' @examples
#' # the number of classes
#' m <- msr_number(landscape, scale = "landscape")
#' msk_global(m)$number
#'
#' # to count the patches of each class, number the patches first; scale =
#' # "class" looks for the table mdf_componentise writes
#' p <- mdf_componentise(landscape, connectivity = 8L, layer = "cover",
#'                       add = "patch")
#' m <- msr_number(p, scale = "class")
#' msk_categories(m)$number
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter
#' @export

msr_number <- function(obj = NULL, scale = "landscape", layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = scale, choices = c("landscape", "class"))
  assertCharacter(x = layer, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  uVals <- sort(unique(vals[!is.na(vals)]))

  if(scale == "landscape"){

    # number of distinct classes
    obj@global[[layer]]$number <- length(uVals)

  } else {

    # number of patches per class, from the patches numbered by
    # mdf_componentise; a class without patches (0) counts none
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
