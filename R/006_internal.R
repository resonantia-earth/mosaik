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

# --- Class tables ------------------------------------------------------------

#' Store a measured value in the class table of a layer
#'
#' The class table lists the classes of a layer under \code{gid}; a value is
#' stored in that order. A layer without class table starts one with the
#' measured classes. A class of the table that was not measured gets
#' \code{NA}; a matrix (one row and column per class) is reordered on both
#' sides.
#'
#' @param obj a mosaik
#' @param layer character(1) the layer
#' @param name character(1) the name of the value
#' @param gid the measured classes
#' @param value one value per class, or a class by class matrix
#' @return the mosaik
#' @noRd

.store_class <- function(obj, layer, name, gid, value) {

  tbl <- obj@categories[[layer]]
  if (is.null(tbl$gid)) tbl <- c(list(gid = gid), tbl)
  idx <- match(tbl$gid, gid)

  if (is.matrix(value)) {
    value <- value[idx, idx, drop = FALSE]
    dimnames(value) <- list(tbl$gid, tbl$gid)
  } else {
    value <- value[idx]
  }
  tbl[[name]] <- value
  obj@categories[[layer]] <- tbl
  obj
}


#' The names of the values measured on a layer
#'
#' Everything in a class table except the class codes, labels and colours,
#' which describe the classes. A table without classes holds fields another
#' package attached, so nothing in it counts.
#' @param entry the class table of one layer.
#' @return character
#' @noRd

.measured_names <- function(entry) {
  if (is.null(entry$gid)) return(character())
  setdiff(names(entry), c("gid", "val", "colour"))
}


#' Warn that a change of the grid leaves measured values behind
#'
#' The values are kept, so they can be compared, but they still describe the
#' map before the change, and a later msr() would read them as current.
#' @param obj the mosaik before the change.
#' @param fn the name of the function.
#' @noRd

.warn_measured <- function(obj, fn) {
  m <- lapply(obj@categories, .measured_names)
  for (nm in names(obj@global)) m[[nm]] <- c(m[[nm]], names(obj@global[[nm]]))
  m <- m[lengths(m) > 0]
  if (!length(m)) return(invisible())
  what <- paste0(names(m), " (", vapply(m, paste, character(1), collapse = ", "),
                 ")", collapse = "; ")
  warning(fn, "() keeps the values measured before it: ", what,
          ". They describe the map before the change; measure again to update ",
          "them.", call. = FALSE)
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

