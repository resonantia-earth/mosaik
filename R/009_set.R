#' Write a layer into a mosaik
#'
#' Return a copy of \code{obj} with one layer's values replaced (or a new layer
#' added), carrying the layer's semantics and the object's provenance forward
#' correctly. This is the single supported way for a function — in this package
#' or a downstream one — to put values into a mosaik.
#'
#' @details
#' \strong{Continuous or categorical.} Passing \code{gid} and \code{val} writes a
#' categorical layer: the category table is set wholesale from them, and only
#' those \code{gid}s that actually occur in \code{values} are registered.
#' Omitting them writes a continuous layer and leaves any existing category
#' table for that layer cleared.
#'
#' \strong{Roles carry forward, unless the layer changes kind.} A value-only
#' operator (\code{\link{mdf_scale}}, \code{\link{mdf_perturb}}, ...) makes no
#' new semantic claim, so it must not silently erase an existing role tag; such
#' callers pass nothing and the prior role is kept. A caller that genuinely
#' changes what the layer *is* — \code{\link{mdf_binarise}} turning a surface
#' into a mask, \code{\link{mdf_componentise}} turning it into patch IDs — passes
#' \code{keep_role = FALSE}. A caller that sets a new meaning passes
#' \code{role}.
#'
#' \strong{\code{tau} rides with the role.} A layer carrying a \code{tau} is one
#' the simulation engine may evolve toward its ambient equilibrium; a layer
#' without one is inert and is never guessed at. It follows the same
#' keep-vs-set discipline as \code{role}, and is dropped when the role is.
#'
#' @param obj [mosaik][mosaik]\cr the mosaik to write into.
#' @param layer [character(1)][character]\cr name of the layer to write. A name
#'   not already present adds a new layer.
#' @param values [numeric(.)][numeric]\cr the new cell values, one per cell.
#'   Compressed with \code{\link[base]{rle}} automatically when that is smaller.
#' @param prov [list][list]\cr a provenance entry from \code{\link{msk_prov}},
#'   appended to \code{obj}'s history.
#' @param role [character(1)][character]\cr role tag to give the layer (e.g.
#'   \code{"surface"}, \code{"temperature"}). \code{NULL} keeps whatever the
#'   layer had, subject to \code{keep_role}.
#' @param keep_role [logical(1)][logical]\cr carry the layer's prior role (and
#'   \code{tau}) forward when overwriting in place. \code{FALSE} for an operator
#'   that changes the layer's kind. Default \code{TRUE}.
#' @param tau [numeric(1)][numeric]\cr relaxation rate: how many timesteps the
#'   engine takes to close the gap between this field and its equilibrium.
#'   \code{NA} clears it explicitly; \code{NULL} keeps the prior value.
#' @param gid [integer(.)][integer]\cr categorical: the group IDs occurring in
#'   \code{values}.
#' @param val [character(.)][character]\cr categorical: the label for each
#'   \code{gid}.
#' @return The mosaik, with the layer written and \code{prov} appended.
#' @seealso \code{\link{msk_prov}} for the provenance entry, \code{\link{mosaik}}
#'   for creating an object in the first place, and \code{\link{msk_add}} for
#'   taking layers from another mosaik.
#' @examples
#' m <- mosaik(extent = c(0, 10, 0, 10), res = 1)
#' v <- runif(100)
#'
#' # a continuous layer, tagged with what it means
#' m <- msk_set(m, "elevation", v * 800, role = "surface",
#'              prov = msk_prov("example", list(range = c(0, 800))))
#'
#' # a categorical one
#' m <- msk_set(m, "cover", rep(1:2, each = 50),
#'              gid = 1:2, val = c("forest", "crop"),
#'              prov = msk_prov("example", list()))
#' @export

msk_set <- function(obj, layer, values, prov = NULL, role = NULL,
                    keep_role = TRUE, tau = NULL, gid = NULL, val = NULL) {

  new_layers <- obj@layers
  # auto-compress
  rv <- rle(values)
  if (utils::object.size(rv) < utils::object.size(values)) {
    new_layers[[layer]] <- rv
  } else {
    new_layers[[layer]] <- values
  }

  categorical <- !is.null(gid)

  # a categorical write sets the whole category entry below, so nothing is
  # carried forward from the layer's previous life
  if (categorical) keep_role <- FALSE

  # carry the prior role forward when overwriting a layer in place: a value-only
  # operator makes no new semantic claim, so it must not silently erase an
  # existing role tag.
  prior_role <- obj@categories[[layer]]$role
  if (is.null(role) && keep_role) role <- prior_role

  # tau rides with the role: a value-only operator keeps the field evolvable, an
  # operator that changes the layer's kind drops it along with the role. NA is
  # the caller explicitly clearing the rate and must not be overridden by the
  # prior value, so only NULL ("caller said nothing") keeps.
  prior_tau <- obj@categories[[layer]]$tau
  if (is.null(tau) && keep_role) tau <- prior_tau
  if (!is.null(tau) && (length(tau) != 1 || is.na(tau))) tau <- NULL

  # clear only the affected layer's categories
  new_categories <- obj@categories
  new_categories[[layer]] <- NULL

  if (categorical) {
    # register only those gids that actually occur in the output
    present <- gid %in% unique(values[!is.na(values)])
    entry <- list(gid = gid[present], val = val[present])
    if (!is.null(role)) entry$role <- role
    if (!is.null(tau))  entry$tau  <- tau
    new_categories[[layer]] <- entry
  } else if (!is.null(role)) {
    # role tag on a continuous layer: no gid/val, just the role
    new_categories[[layer]] <- list(role = role)
    if (!is.null(tau)) new_categories[[layer]]$tau <- tau
  }

  new_mosaik(
    extent     = obj@extent,
    dims       = obj@dims,
    layers     = new_layers,
    categories = new_categories,
    patches    = obj@patches,
    global     = obj@global,
    crs        = obj@crs,
    provenance = if (is.null(prov)) obj@provenance else c(obj@provenance, list(prov))
  )
}
