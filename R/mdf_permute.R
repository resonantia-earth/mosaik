#' Permute the values in a mosaik
#'
#' The permutation of a set of cell values leads to a systematic
#' (re)arrangement of all the members of the set.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param type [character(1)][character]\cr the permutation type. Either
#'   \code{"invert"}, \code{"revert"}, \code{"descending"}, \code{"ascending"}
#'   or \code{"cycle"}.
#' @param by [integerish(1)][integer]\cr value by which to apply the
#'   permutation; only for \code{type = "cycle"}.
#' @param layer [character(1)][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}.
#' @examples
#' mdf_permute(landscape, type = "invert", layer = "cover")
#' mdf_permute(landscape, type = "revert", layer = "cover")
#' mdf_permute(landscape, type = "ascending", layer = "intensity")
#' mdf_permute(landscape, type = "descending", layer = "intensity")
#' mdf_permute(landscape, type = "cycle", by = 2L, layer = "cover")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertChoice assertIntegerish
#'   assertCharacter
#' @export

mdf_permute <- function(obj = NULL,
                        type = "invert",
                        by = NULL,
                        layer = NULL,
                        add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertChoice(x = type, choices = c("invert", "revert", "descending", "ascending", "cycle"))
  if(type == "cycle"){
    assertIntegerish(x = by)
  }
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  uVals <- sort(unique(vals[!is.na(vals)]))

  # body ----
  if(type == "invert"){
    newVals <- max(uVals) - uVals
  } else if(type == "revert"){
    newVals <- rev(uVals)
  } else if(type == "descending"){
    newVals <- sort(uVals, decreasing = TRUE)
  } else if(type == "ascending"){
    newVals <- sort(uVals, decreasing = FALSE)
  } else if(type == "cycle"){
    newVals <- uVals + by
    newVals[newVals > max(uVals)] <- newVals[newVals > max(uVals)] - max(uVals)
    newVals[newVals < min(uVals)] <- newVals[newVals < min(uVals)] + max(uVals)
  }

  temp <- vals
  for(i in seq_along(uVals)){
    temp[vals == uVals[i]] <- newVals[i]
  }

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_permute", list(type = type, by = by, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
