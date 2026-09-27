#' Add random noise to cell values
#'
#' Perturb cell values by adding random noise drawn from a normal distribution
#' with mean 0 and standard deviation \code{sd}.
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param sd [numeric(1)][numeric]\cr standard deviation of the noise.
#'   Defaults to 1.
#' @param layer [character(1)][character]\cr the layer to perturb. Defaults to
#'   the first layer.
#' @param add [character(1)][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions with perturbed cell values. NA cells
#'   remain NA.
#' @examples
#' mdf_perturb(landscape, sd = 10, layer = "intensity")
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertNumber assertCharacter
#' @importFrom stats rnorm
#' @export

mdf_perturb <- function(obj = NULL,
                        sd = 1,
                        layer = NULL,
                        add = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumber(x = sd, lower = 0, finite = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  noise <- rnorm(length(vals), mean = 0, sd = sd)
  temp <- vals + noise
  temp[is.na(vals)] <- NA

  # build output ----
  out_layer <- .resolve_add(obj, layer, add)
  prov <- msk_prov("mdf_perturb", list(sd = sd, layer = out_layer))
  msk_set(obj, out_layer, temp, prov)
}
