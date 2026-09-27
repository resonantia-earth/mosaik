#' Mosaik constructor and validator
#'
#' Low-level constructor for the \code{mosaik} class. Validates the raw slots and
#' returns a valid mosaik; every mosaik-producing function funnels through here.
#' Not user-facing — the user-facing creator is \code{mosaik()}, which dispatches
#' by input type and computes defaults before calling this.
#'
#' @param extent [numeric(4)][numeric]\cr spatial extent as
#'   \code{c(xmin, xmax, ymin, ymax)}.
#' @param dims [integer(2)][integer]\cr grid dimensions as
#'   \code{c(ncols, nrows)}.
#' @param layers [list][list]\cr named list of layer values (vectors or
#'   \code{rle}). If empty, a single \code{NA} layer named \code{"values"} is
#'   created.
#' @param categories [list][list]\cr per-layer category tables. Each element
#'   is a list with at least \code{gid}.
#' @param patches [list][list]\cr patch-level attributes.
#' @param global [list][list]\cr landscape-level metrics.
#' @param crs [character(1)][character]\cr CRS in proj4 notation, or
#'   \code{NA_character_} for cartesian.
#' @param provenance [list][list]\cr provenance history entries.
#' @return An object of class \code{mosaik}.
#' @importFrom checkmate assertNumeric assertIntegerish assertCharacter assertList
#' @importFrom methods new
#' @importFrom utils object.size
#' @noRd

new_mosaik <- function(extent,
                       dims,
                       layers = list(),
                       categories = list(),
                       patches = list(),
                       global = list(),
                       crs = NA_character_,
                       provenance = list()) {

  # --- extent ---
  assertNumeric(x = extent, len = 4, any.missing = FALSE)
  if (extent[2] <= extent[1]) stop("extent xmax must be > xmin.")
  if (extent[4] <= extent[3]) stop("extent ymax must be > ymin.")

  # --- dims ---
  assertIntegerish(x = dims, len = 2, lower = 1, any.missing = FALSE)
  dims <- as.integer(dims)
  n_cells <- prod(dims)

  # --- layers ---
  assertList(x = layers)
  layerNames <- names(layers)
  if (length(layers) > 0 && (is.null(layerNames) || any(layerNames == ""))) {
    stop("'layers' must be a named list.")
  }

  # validate and auto-compress
  for (nm in layerNames) {
    v <- layers[[nm]]
    actual_len <- if (inherits(v, "rle")) sum(v$lengths) else length(v)
    if (actual_len != n_cells) {
      stop(sprintf("Layer '%s' has %d values but grid has %d cells.",
                   nm, actual_len, n_cells))
    }
    if (!inherits(v, "rle")) {
      rv <- rle(v)
      if (object.size(rv) < object.size(v)) {
        layers[[nm]] <- rv
      }
    }
  }

  # --- categories ---
  assertList(x = categories)
  if (length(categories) > 0) {
    for (gnm in names(categories)) {
      if (!gnm %in% layerNames) {
        stop("categories entry '", gnm, "' does not match any layer. ",
             "Available layers: ", paste(layerNames, collapse = ", "), ".")
      }
      entry <- categories[[gnm]]
      if (!is.list(entry) || is.data.frame(entry)) {
        stop("Categories entry '", gnm, "' must be a list (e.g. list(gid = ..., val = ...)).")
      }
      # role-only entries (continuous layers tagged with a role but no
      # categorical values) are valid — skip gid cross-check for those
      if (!"gid" %in% names(entry)) {
        if ("role" %in% names(entry)) next
        stop("Categories entry '", gnm, "' must contain 'gid'.")
      }
      # cross-check gid vs cell values
      layerVals <- layers[[gnm]]
      if (inherits(layerVals, "rle")) layerVals <- inverse.rle(layerVals)
      uVals <- sort(unique(layerVals[!is.na(layerVals)]))
      missing <- setdiff(uVals, entry$gid)
      if (length(missing) > 0) {
        warning("Layer '", gnm, "' contains cell values (",
                paste(missing, collapse = ", "), ") not listed in categories$gid.")
      }
      extra <- setdiff(entry$gid, uVals)
      if (length(extra) > 0) {
        warning("Categories entry '", gnm, "' has gid values (",
                paste(extra, collapse = ", "), ") not found in layer cell values.")
      }
    }
  }

  # --- crs ---
  assertCharacter(x = crs, len = 1)

  # --- provenance ---
  if (!is.list(provenance)) provenance <- list(provenance)

  # --- assemble ---
  new(Class      = "mosaik",
      extent     = extent,
      dims       = dims,
      layers     = layers,
      categories = categories,
      patches    = patches,
      global     = global,
      crs        = crs,
      provenance = provenance)
}


