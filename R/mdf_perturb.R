#' Add random noise to cell values
#'
#' Perturb cell values by adding random noise drawn from a normal distribution
#' with mean 0 and standard deviation \code{sd}.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param sd [`numeric(1)`][numeric]\cr standard deviation of the noise.
#'   Defaults to 1.
#' @param layer [`character(1)`][character]\cr the layer to perturb. Defaults to
#'   the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions with perturbed cell values. NA cells
#'   remain NA.
#' @examples
#' # the canopy height with normal noise of standard deviation 2 added
#' m <- mdf_perturb(landscape, sd = 2, layer = "canopy", add = "perturbed")
#' msk_vis(m, .layer("canopy"), .layer("perturbed"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertNumber assertCharacter
#' @importFrom stats rnorm
#' @export

mdf_perturb <- function(obj = NULL,
                        sd = 1,
                        layer = NULL,
                        add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

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
  .update_mosaik(obj, values = temp, step = step)
}
