#' Landscape mosaic class (S4)
#'
#' A mosaik represents a landscape as one or more co-registered grid layers
#' sharing the same extent, resolution, and coordinate reference system. It
#' stores gridded landscape data as a transparent, inspectable S4 object.
#'
#' @slot extent [`numeric(4)`][numeric]\cr \code{c(xmin, xmax, ymin, ymax)}.
#' @slot dims [`integer(2)`][integer]\cr \code{c(ncols, nrows)}.
#' @slot layers [`list`][list]\cr named list of cell value vectors. It stores the
#'   raw cell values. Each layer is a flat vector in row-major order
#'   (RLE-compressed when beneficial).
#' @slot categories [`list`][list]\cr named list of per-layer category tables. It
#'   provides metadata for categorical layers. Each entry is a list keyed by
#'   layer name, containing at least \code{gid} (the class IDs present in the
#'   layer). Optional elements include \code{val} (human-readable labels such as
#'   \code{"forest"} or \code{"water"}) and class-level metrics added by
#'   functions such as \code{\link{msr_area}(scale = "class")}.
#' @slot patches [`list`][list]\cr named list of per-layer patch tables. An entry
#'   is created by \code{\link{mdf_componentise}}: the layer holding the patch
#'   numbers (\code{ids}), the \code{connectivity}, and the \code{class} and
#'   number (\code{patch}) of each patch. Patch-level measures such as
#'   \code{\link{msr_area}(scale = "patch")} add their values to it.
#' @slot global [`list`][list]\cr named list of per-layer landscape-level values,
#'   from functions such as \code{\link{msr_area}(scale = "landscape")}.
#' @slot crs [`character(1)`][character]\cr coordinate reference system string or
#'   \code{NA_character_}.
#' @slot provenance [`list`][list]\cr processing history, i.e. every operation
#'   applied to the mosaik, enabling full traceability of the processing chain.
#' @exportClass mosaik

setClass(Class = "mosaik",
         slots = c(extent     = "numeric",
                   dims       = "integer",
                   layers     = "list",
                   patches    = "list",
                   categories = "list",
                   global     = "list",
                   crs        = "character",
                   provenance = "list"
         )
)

setValidity("mosaik", function(object){

  errors <- character()

  # extent
  if (length(object@extent) != 4) {
    errors <- c(errors, "'extent' must be length 4 (xmin, xmax, ymin, ymax).")
  } else {
    if (object@extent[2] <= object@extent[1]) {
      errors <- c(errors, "extent xmax must be greater than xmin.")
    }
    if (object@extent[4] <= object@extent[3]) {
      errors <- c(errors, "extent ymax must be greater than ymin.")
    }
  }

  # dims
  if (length(object@dims) != 2) {
    errors <- c(errors, "'dims' must be length 2 (ncols, nrows).")
  } else {
    if (any(object@dims < 1L)) {
      errors <- c(errors, "'dims' must be positive integers.")
    }
  }

  # layers
  if (!is.list(object@layers)) {
    errors <- c(errors, "'layers' must be a named list.")
  } else if (length(object@layers) > 0 &&
             (is.null(names(object@layers)) || any(names(object@layers) == ""))) {
    errors <- c(errors, "'layers' must be a named list.")
  } else {
    n_cells <- prod(object@dims)
    for (nm in names(object@layers)) {
      v <- object@layers[[nm]]
      actual_len <- if (inherits(v, "rle")) sum(v$lengths) else length(v)
      if (actual_len != n_cells) {
        errors <- c(errors, sprintf("layer '%s' has %d values but grid has %d cells.",
                                    nm, actual_len, n_cells))
      }
    }
  }

  # patches and global: one list of results per layer, like categories
  for (s in c("patches", "global")) {
    tbl <- methods::slot(object, s)
    if (!is.list(tbl)) {
      errors <- c(errors, sprintf("'%s' must be a list.", s))
      next
    }
    for (nm in names(tbl)) {
      if (!nm %in% names(object@layers)) {
        errors <- c(errors, sprintf("%s entry '%s' does not match any layer.", s, nm))
      } else if (!is.list(tbl[[nm]]) || is.data.frame(tbl[[nm]])) {
        errors <- c(errors, sprintf("%s entry '%s' must be a list.", s, nm))
      }
    }
  }

  # categories
  if (!is.list(object@categories)) {
    errors <- c(errors, "'categories' must be a list.")
  } else {
    for (nm in names(object@categories)) {
      entry <- object@categories[[nm]]
      if (!nm %in% names(object@layers)) {
        errors <- c(errors, sprintf("categories entry '%s' does not match any layer.", nm))
      } else if (!is.list(entry) || is.data.frame(entry)) {
        errors <- c(errors, sprintf("categories entry '%s' must be a list.", nm))
      } else if ("val" %in% names(entry) && !"gid" %in% names(entry)) {
        # an entry without gid belongs to a continuous layer and carries fields
        # another package attached to it; only labels without IDs are an error
        errors <- c(errors, sprintf("categories entry '%s' has 'val' but no 'gid'.", nm))
      }
    }
  }

  # crs
  if (!is.character(object@crs) || length(object@crs) != 1) {
    errors <- c(errors, "'crs' must be a single character string.")
  }

  # provenance
  if (!is.list(object@provenance)) {
    errors <- c(errors, "'provenance' must be a list.")
  }

  if (length(errors) == 0) TRUE else errors
})
