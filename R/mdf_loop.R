#' Repeat a recipe until a condition is met
#'
#' Apply a recorded sequence of steps (a recipe, see \code{\link{mdf}}) over and
#' over, either a fixed number of times or until a condition holds. Iteration is
#' a separate operator rather than an argument of the operators that happen to
#' need it, so every \code{mdf_*} keeps doing one thing.
#'
#' @param obj [mosaik]\cr the mosaik to modify.
#' @param recipe [mosaik]\cr a recipe shell holding the steps of one iteration,
#'   built by calling \code{mdf_*} with \code{obj = NULL}.
#' @param times [integerish(1)][integer]\cr the maximum number of iterations.
#'   \code{Inf} (default) runs until a stopping condition is met, and then one
#'   of \code{until} or \code{stable} is required.
#' @param until an expression over the layers, evaluated after each iteration;
#'   iteration stops once it is \code{TRUE}. It must reduce to a single value,
#'   e.g. \code{sum(corefill) == 0}.
#' @param stable [logical(1)][logical]\cr stop once an iteration no longer
#'   changes the layer named by \code{layer} -- a fixpoint. This is the
#'   condition that cannot be written as an expression over layers, because it
#'   compares against the previous state.
#' @param layer [character(1)][character]\cr the layer watched by
#'   \code{stable}. Defaults to the first layer.
#' @return A mosaik of the same dimensions as \code{obj}.
#' @details
#'   The three stopping conditions cover the three shapes of loop that occur in
#'   raster algorithms:
#'
#'   \describe{
#'     \item{\code{times = n}}{a for-loop: repeat exactly \code{n} times, e.g. a
#'       geodesic dilation of a given size.}
#'     \item{\code{until = <expr>}}{a while-loop with an explicit condition,
#'       e.g. peel a set until nothing is left.}
#'     \item{\code{stable = TRUE}}{run to a fixpoint, e.g. spread a marker
#'       through a mask until it stops growing (reconstruction by dilation).}
#'   }
#'
#'   \code{times} may be combined with either condition, in which case it caps
#'   the iterations and guards against a loop that never settles.
#' @examples
#' # the forest (class 47), and one of its cells as a seed
#' m <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest")
#' seed <- as.numeric(seq_len(msk_ncells(m)) ==
#'                    which(msk_pull(m, "forest") == 1)[1])
#' m <- msk_set(m, "seed", seed, prov = msk_prov("seed", list()))
#'
#' # spread the seed through the forest until it stops growing: what is left is
#' # the forest patch the seed sits in
#' grow <- mdf_dilate(struct = msk_struct("square", width = 3), layer = "seed") |>
#'   mdf_filter(seed == 1 & forest == 1, value = TRUE, background = 0,
#'              add = "seed")
#' m <- mdf_loop(m, grow, stable = TRUE, layer = "seed")
#' sum(msk_pull(m, "seed"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertNumber assertFlag assertCharacter
#' @export

mdf_loop <- function(obj = NULL,
                     recipe,
                     times = Inf,
                     until = NULL,
                     stable = FALSE,
                     layer = NULL){

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertClass(x = recipe, classes = "mosaik")
  assertNumber(x = times, lower = 1)
  assertFlag(x = stable)
  assertCharacter(x = layer, len = 1, null.ok = TRUE)

  condition <- substitute(until)
  if (is.name(condition) &&
      exists(as.character(condition), envir = parent.frame(), inherits = FALSE)) {
    bound <- get(as.character(condition), envir = parent.frame(), inherits = FALSE)
    if (is.call(bound) || is.name(bound)) condition <- bound
  }
  hasUntil <- !is.null(condition) && !identical(condition, quote(NULL))

  if (is.infinite(times) && !hasUntil && !stable) {
    stop("'times = Inf' needs a stopping condition: pass 'until' or ",
         "'stable = TRUE'.", call. = FALSE)
  }

  if (is.null(layer)) layer <- names(obj@layers)[1]

  # body ----
  i <- 0
  while (i < times) {
    before <- if (stable) msk_pull(obj, layer) else NULL

    obj <- mdf(obj, recipe)
    i <- i + 1

    if (stable) {
      after <- msk_pull(obj, layer)
      if (isTRUE(all.equal(before, after))) break
    }

    if (hasUntil) {
      env <- lapply(names(obj@layers), function(nm) msk_pull(obj, nm))
      names(env) <- names(obj@layers)
      done <- eval(condition, envir = env, enclos = parent.frame())
      if (length(done) != 1 || !is.logical(done)) {
        stop("'until' must reduce to a single logical value.", call. = FALSE)
      }
      if (isTRUE(done)) break
    }
  }

  # build output ----
  # .make_prov already returns a named one-element list; wrapping it again
  # would drop the name
  prov <- msk_prov("mdf_loop", list(iterations = i, layer = layer))
  obj@provenance <- c(obj@provenance, prov)
  obj
}
