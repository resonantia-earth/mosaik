# Internal helpers used by mdf_*, msr_* and msk_* functions.
# Not exported.

#' An empty, valueless mosaik used to record a recipe
#'
#' Geometry is a placeholder (1x1) and is ignored on replay; only
#' \code{@provenance} carries the recorded calls.
#' @return a 1x1 mosaik with no layers.
#' @noRd
.recipe_shell <- function() {
  new_mosaik(extent = c(0, 1, 0, 1), dims = c(1L, 1L))
}

#' Test whether an object is a recipe shell
#'
#' TRUE if \code{obj} is a recipe shell (no layers), i.e. a call made to record
#' rather than execute.
#' @param obj a mosaik, or NULL.
#' @return logical(1).
#' @noRd
.is_recipe <- function(obj) {
  is.null(obj) || (methods::is(obj, "mosaik") && length(obj@layers) == 0)
}

#' Record one step of a recipe
#'
#' Append the captured call (function name plus the arguments other than
#' \code{obj}) to the shell's provenance and return the shell. The call is
#' captured by the calling \code{mdf_*} via \code{match.call()}.
#' @param obj a recipe shell, or NULL (a fresh shell is created).
#' @param cl the captured call (from \code{match.call()}).
#' @return the recipe shell with one step appended to \code{@provenance}.
#' @details Arguments are evaluated as they are recorded, so that a step stores
#'   the object a call produced (\code{struct = msk_struct(...)}) rather than
#'   the call itself, and replaying with \code{do.call} needs no environment.
#'   Predicate arguments are the exception: they name layers of the mosaik the
#'   recipe will later be applied to, so they cannot be evaluated now and are
#'   stored unevaluated. \code{.quoted_args} lists them per function.
#' @noRd

# arguments that must NOT be evaluated when a step is recorded, because they
# are expressions over the layers of the future target mosaik
.quoted_args <- list(mdf_filter = "expr",
                     mdf_loop   = "until")

.record_step <- function(obj, cl) {
  if (is.null(obj)) obj <- .recipe_shell()
  fn <- as.character(cl[[1]])
  args <- as.list(cl)[-1]
  args[["obj"]] <- NULL
  keep <- .quoted_args[[fn]]
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  for (i in seq_along(args)) {
    if (nms[i] %in% keep) next
    a <- args[[i]]
    if (is.name(a) || is.call(a)) args[[i]] <- eval(a, parent.frame(2))
  }
  obj@provenance <- c(obj@provenance, msk_prov(fn, args, step = TRUE))
  obj
}

#' Resolve a layer-name argument to a values vector
#'
#' A second layer is always referred to **by name**, never by passing another
#' mosaik: a layer from elsewhere is brought in with \code{msk_add()} first, so
#' the grid-conformity check lives in one place and recipe steps can name the
#' layer they depend on. Used by \code{mdf_mask} / \code{mdf_zonal} /
#' \code{mdf_layerise}.
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

# --- Patch decomposition cache --------------------------------------------

#' Ensure a global patch-ID layer exists
#'
#' Checks whether an internal layer \code{"_patches"} already exists (and was
#' computed from the same cover \code{layer}). If not, loops over unique class
#' values, binarises and calls \code{componentsCpp} for each class, and assigns
#' globally unique patch IDs: class 1's patches get IDs 1..n1, class 2's
#' patches get (n1+1)..(n1+n2), etc. The result is stored as a single integer
#' layer \code{"_patches"}.
#'
#' Combined with the original cover layer (which gives each cell's class), any
#' downstream function can reconstruct per-class patch membership trivially.
#'
#' Also populates \code{@patches$class} and \code{@patches$patch} (the flat
#' roster) as a side-effect.
#'
#' @param obj a mosaik
#' @param layer character(1) the cover layer to decompose
#' @return the mosaik, potentially with \code{"_patches"} added
#' @noRd

.ensure_patch_layer <- function(obj, layer) {

  # check cache: reuse if present and from the same layer
  if ("_patches" %in% names(obj@layers) &&
      identical(obj@patches$._patch_layer, layer)) {
    return(obj)
  }

  vals <- msk_pull(obj, layer)
  dims <- obj@dims
  uVals <- sort(unique(vals[!is.na(vals)]))

  patch_layer <- rep(NA_integer_, length(vals))
  offset <- 0L
  all_classes <- NULL
  all_patches <- NULL

  for (i in seq_along(uVals)) {
    temp_vals <- vals
    temp_vals[temp_vals != uVals[i]] <- NA
    temp_cc <- componentsCpp(vals = temp_vals, nrow = dims[2], ncol = dims[1])

    local_ids <- sort(unique(temp_cc[!is.na(temp_cc)]))
    n_local <- length(local_ids)

    # write globally offset IDs into the patch layer
    cells <- which(!is.na(temp_cc))
    patch_layer[cells] <- temp_cc[cells] + offset

    all_classes <- c(all_classes, rep(uVals[i], n_local))
    all_patches <- c(all_patches, local_ids + offset)

    offset <- offset + n_local
  }

  # store the layer
  new_layers <- obj@layers
  new_layers[["_patches"]] <- patch_layer
  obj <- new_mosaik(
    extent     = obj@extent,
    dims       = obj@dims,
    layers     = new_layers,
    categories = obj@categories,
    patches    = obj@patches,
    global     = obj@global,
    crs        = obj@crs,
    provenance = obj@provenance
  )

  # store the flat roster and cache key
  obj@patches$class <- all_classes
  obj@patches$patch <- all_patches
  obj@patches$._patch_layer <- layer

  obj
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

.resolve_add <- function(obj, layer, add) {
  if (is.null(add)) {
    if (is.null(layer) || length(layer) == 0 || is.na(layer)) return("values")
    return(layer)
  }
  add
}

