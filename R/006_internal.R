# Internal helpers used by mdf_*, msr_* and msk_* functions.
# Not exported.

#' An empty, valueless mosaik used to record a recipe
#'
#' Geometry is a placeholder (1x1) and is ignored on replay; only
#' \code{@provenance} carries the recorded calls.
#' @return a 1x1 mosaik with no layers.
#' @noRd
.recipe_shell <- function() {
  methods::new("mosaik", extent = c(0, 1, 0, 1), dims = c(1L, 1L),
               crs = NA_character_)
}

#' Test whether an object is a recipe shell
#'
#' A recipe shell is a mosaik without layers: a call on it, or on \code{NULL},
#' records the step instead of running it.
#' @param obj a mosaik, or NULL.
#' @return logical(1).
#' @noRd
.is_recipe_shell <- function(x) {
  methods::is(x, "mosaik") && length(x@layers) == 0
}

.is_recipe <- function(obj) is.null(obj) || .is_recipe_shell(obj)

#' Resolve a layer-name argument to a values vector
#'
#' A second layer is always referred to **by name**, never by passing another
#' mosaik: a layer from elsewhere is brought in with \code{msk_add()} first, so
#' the grid-conformity check lives in one place and recipe steps can name the
#' layer they depend on. Used by \code{mdf_summarise}.
#' @param obj the mosaik holding the layer.
#' @param x a character layer name in \code{obj}.
#' @param arg the argument name, for the error message.
#' @return a values vector.
#' @noRd
.pull_layer <- function(obj, x, arg = "by") {
  if (!is.character(x) || length(x) != 1) {
    stop("'", arg, "' must be the name of a layer in 'obj'. To use a layer of ",
         "another mosaik, add it with msk_add() first.", call. = FALSE)
  }
  if (!(x %in% names(obj@layers))) {
    stop("layer '", x, "' not found in 'obj'.", call. = FALSE)
  }
  msk_pull(obj, x)
}

#' Resolve a summary shorthand to a function
#'
#' The one vocabulary of summaries shared by \code{mdf_summarise} (over the
#' cells of a zone) and \code{mdf_blend} (over the layers of a cell).
#' @param fun a shorthand or a function.
#' @param caller the calling function, for the error message.
#' @return a function taking a numeric vector and returning one value.
#' @noRd
.summary_fun <- function(fun, caller) {
  if (is.function(fun)) return(fun)
  if (!is.character(fun) || length(fun) != 1) {
    stop("'fun' must be a shorthand or a function; see ?", caller, ".",
         call. = FALSE)
  }
  switch(fun,
         # the extreme of no values is unknown, not -Inf or Inf
         "max"        = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE),
         "min"        = function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
         "sum"        = function(x) sum(x, na.rm = TRUE),
         "mean"       = function(x) mean(x, na.rm = TRUE),
         "median"     = function(x) stats::median(x, na.rm = TRUE),
         "any"        = function(x) as.numeric(any(x != 0, na.rm = TRUE)),
         "all"        = function(x) as.numeric(all(x != 0, na.rm = TRUE)),
         "n"          = function(x) length(x),
         "n_distinct" = function(x) length(unique(x[!is.na(x) & x != 0])),
         "unique"     = function(x) {
           u <- unique(x[!is.na(x) & x != 0])
           if (length(u) == 1) u else NA_real_
         },
         stop("unknown 'fun' shorthand '", fun, "'; see ?", caller, ".",
              call. = FALSE))
}

#' Trace a path from target to source using predecessor map
#'
#' Walks backward through the predecessor vector returned by
#' \code{costDistanceCpp} to reconstruct the shortest path.
#'
#' @param pred integer vector of 1-based predecessor indices (-1 = source,
#'   0 = unreachable).
#' @param target integer(1), 1-based target cell index.
#' @return Integer vector of 1-based cell indices along the path (from target
#'   to source, inclusive).
#' @noRd

.trace_path <- function(pred, target) {
  path <- integer(0)
  cur <- target

  max_steps <- length(pred)  # safety limit
  steps <- 0L

  while (TRUE) {
    steps <- steps + 1L
    if (steps > max_steps) {
      warning("Path tracing exceeded maximum steps, possible cycle.")
      break
    }

    path <- c(path, cur)

    p <- pred[cur]
    if (p == -1L) break    # reached source
    if (p == 0L) {
      warning(sprintf("Target cell %d is unreachable from source.", target))
      break
    }

    cur <- p
  }

  rev(path)  # return source-to-target order
}

# --- Patches ---------------------------------------------------------------

#' The patches of a layer, as numbered by mdf_componentise
#'
#' Patches are never found by a measure: \code{\link{mdf_componentise}} numbers
#' them, with a connectivity the user states, and records under
#' \code{@patches[[layer]]} the layer holding the numbers (\code{ids}) and the
#' class and number of each patch. The patch-level
#' measures read that record here, and stop with a pointer to
#' \code{mdf_componentise} if there is none.
#'
#' @param obj a mosaik
#' @param layer character(1) the layer whose patches are measured
#' @return list: \code{ids} (patch number per cell), \code{class} and
#'   \code{patch} (one entry per patch, in record order).
#' @noRd

.patches_of <- function(obj, layer) {

  rec <- obj@patches[[layer]]
  if (is.null(rec$ids)) {
    stop(sprintf(paste0(
      "Layer '%s' has no patches. Number them first with ",
      "mdf_componentise(layer = \"%s\", ...), which sets how cells connect ",
      "into patches."), layer, layer), call. = FALSE)
  }
  if (!rec$ids %in% names(obj@layers)) {
    stop(sprintf(paste0(
      "The patches of layer '%s' were numbered in layer '%s', which is no ",
      "longer in 'obj'. Number them again with mdf_componentise()."),
      layer, rec$ids), call. = FALSE)
  }
  ids <- msk_pull(obj, rec$ids)

  # patches the map border cuts: their full extent is unknown
  grid <- matrix(ids, nrow = obj@dims[2], ncol = obj@dims[1], byrow = TRUE)
  border <- c(grid[1, ], grid[nrow(grid), ], grid[, 1], grid[, ncol(grid)])

  list(ids = ids, class = rec$class, patch = rec$patch,
       clipped = rec$patch %in% border)
}


#' Unlink patch records from layers that are no longer present
#'
#' After layers are selected or removed, a patch record whose number layer
#' went with them keeps its measured values but loses the link, so a further
#' patch-level measure asks for \code{mdf_componentise} again.
#'
#' @param patches the \code{@patches} slot, already subset to the kept layers.
#' @param layers the names of the layers that remain.
#' @return the patches slot.
#' @noRd

.unlink_patches <- function(patches, layers) {
  for (nm in names(patches)) {
    if (!is.null(patches[[nm]]$ids) && !patches[[nm]]$ids %in% layers) {
      patches[[nm]]$ids <- NULL
    }
  }
  patches
}


#' Resolve layer name for mdf_* functions
#'
#' When \code{add = NULL}, returns \code{layer} (overwrite mode).
#' When \code{add} is a character string, returns that string as the new layer
#' name (accumulation mode).
#' @param obj mosaik object
#' @param layer resolved input layer name
#' @param add NULL (overwrite) or character(1) (new layer name)
#' @return character(1) layer name
#' @noRd

