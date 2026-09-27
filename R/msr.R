#' Measure a composite metric
#'
#' Compute a derived landscape metric from an equation over existing
#' primitives measured by \code{msr_*} functions. Variables in the equation use
#' \code{metric.scale} notation to say which attribute table each value comes
#' from; the result is stored automatically at the level implied by its length.
#'
#' Provenance is recorded with the equation, so that \code{sim_target()}
#' can reference derived metrics by name and \code{sim()} can look up
#' the equation for decomposition.
#'
#' @param obj [mosaik]\cr the mosaik with primitives already computed.
#' @param equation [character(1)][character]\cr a mathematical expression using
#'   \code{metric.scale} notation, where scale is one of \code{class},
#'   \code{patch}, \code{landscape}, or \code{cell}. See Details for the
#'   \code{cell} scale, and the examples for equation forms.
#' @param label [character(1)][character]\cr name for the derived metric.
#'   Becomes the column name in the attribute table.
#' @param layer [character(1)][character]\cr which layer's attribute table
#'   to use. Defaults to the first layer.
#' @return The mosaik with the derived metric added to the attribute table
#'   and provenance recorded.
#' @details
#'   The \code{cell} scale reads a per-cell internal layer named
#'   \code{"_<metric>"} (e.g. \code{distance.cell} reads \code{"_distance"},
#'   written by \code{msr_cost(scale = "cell")}), groups its values by patch,
#'   and reduces them to a scalar per patch.
#'
#'   Most equations combine scalar or per-class primitives and evaluate directly.
#'   One case is special: a \strong{patch-distance} primitive
#'   (\code{distance.patch}) is not a single value per patch but a whole
#'   \emph{matrix} — for every patch, its distance to every other patch of the
#'   same class. When the equation uses it, \code{msr()} evaluates the equation
#'   once per patch: \code{distance.patch} becomes that patch's row of the matrix
#'   (its distances to all the others), while any other \code{.patch} variable
#'   becomes the matching value for that same patch. A patch's distance to itself
#'   is \code{Inf}, so self-terms drop out on their own (e.g. \code{1/Inf = 0}).
#'   Each per-patch equation must return a single number — that is what gets
#'   stored. So \code{"min(distance.patch)"} yields each patch's
#'   nearest-neighbour distance.
#' @family measure
#' @examples
#' m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
#'             vals = list(cover = sample(1:3, 100, replace = TRUE)))
#' m <- msr_area(m, scale = "class")
#' m <- msr_perimeter(m, scale = "class")
#'
#' # a ratio of two per-class primitives
#' m <- msr(m, equation = "perimeter.class / area.class",
#'          label = "edge_density")
#'
#' # mixing a class primitive with a landscape one
#' m <- msr_area(m, scale = "landscape")
#' m <- msr(m, equation = "area.class / area.landscape * 100",
#'          label = "prop_area")
#'
#' # a per-patch reduction of the patch-distance matrix. A patch-scale
#' # primitive has to come first: it is what decomposes the layer into patches
#' # and creates the table the result is stored in.
#' m <- msr_area(m, scale = "patch")
#' m <- msr_cost(m, scale = "patch")
#' m <- msr(m, equation = "min(distance.patch)", label = "enn")
#' @importFrom checkmate assertCharacter assertClass
#' @export

msr <- function(obj, equation, label, layer = NULL) {

  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = equation, len = 1)
  assertCharacter(x = label, len = 1)
  assertCharacter(x = layer, null.ok = TRUE)

  if (is.null(layer)) layer <- names(obj@layers)[1]

  scale_map <- list(
    class     = obj@categories[[layer]],
    patch     = obj@patches,
    landscape = obj@global
  )

  # parse equation variables
  eq_parsed <- parse(text = equation)
  eq_vars <- all.vars(eq_parsed)

  # resolve each variable via metric.scale notation
  env <- list()
  scales_used <- character(0)
  has_matrix_patch <- FALSE
  matrix_patch_vars <- character(0)
  has_cell <- FALSE
  cell_vars <- character(0)

  for (v in eq_vars) {
    parts <- strsplit(v, ".", fixed = TRUE)[[1]]

    if (length(parts) != 2 || !(parts[2] %in% c("class", "patch", "landscape", "cell"))) {
      stop(sprintf(
        "Variable '%s' must use metric.scale notation (e.g. 'area.class', 'distance.cell').",
        v))
    }

    metric <- parts[1]
    v_scale <- parts[2]

    if (v_scale == "cell") {
      # cell scale: look up internal layer named "_<metric>"
      cell_layer <- paste0("_", metric)
      if (!(cell_layer %in% names(obj@layers))) {
        stop(sprintf("Internal layer '%s' not found. Run a measure at scale = 'cell' that writes it first.",
                     cell_layer))
      }
      env[[v]] <- msk_pull(obj, cell_layer)
      has_cell <- TRUE
      cell_vars <- c(cell_vars, v)
    } else {
      tbl <- scale_map[[v_scale]]

      if (is.null(tbl) || is.null(tbl[[metric]])) {
        stop(sprintf("Metric '%s' not found at %s level for layer '%s'.",
                     metric, v_scale, layer))
      }

      val <- tbl[[metric]]

      # detect matrix-valued patch data (e.g. distance stored as list of matrices)
      if (v_scale == "patch" && is.list(val) && !is.data.frame(val)) {
        has_matrix_patch <- TRUE
        matrix_patch_vars <- c(matrix_patch_vars, v)
      }

      env[[v]] <- val
    }
    scales_used <- c(scales_used, v_scale)
  }

  if (has_cell) {
    result <- .derive_cellwise(obj, env, eq_parsed, cell_vars, scale_map, layer)
  } else if (has_matrix_patch) {
    result <- .derive_rowwise(env, eq_parsed, matrix_patch_vars, scale_map)
  } else {
    result <- eval(eq_parsed, envir = env)
  }

  # determine storage level from result length
  n_class <- length(scale_map$class$gid)
  n_patch <- length(scale_map$patch$patch)

  if (has_cell) {
    # cell-wise evaluation always produces per-patch results
    obj@patches[[label]] <- result
    store_scale <- "patch"
  } else if (length(result) == n_class && n_class > 0) {
    obj@categories[[layer]][[label]] <- result
    store_scale <- "class"
  } else if (length(result) == n_patch && n_patch > 0) {
    obj@patches[[label]] <- result
    store_scale <- "patch"
  } else if (length(result) == 1) {
    obj@global[[label]] <- result
    store_scale <- "landscape"
  } else if (n_patch == 0 && any(grepl("\\.patch$", eq_vars))) {
    # the equation asks for patch scale but nothing has decomposed the layer
    # into patches yet, so there is no table to store the result in. This is a
    # missing step rather than a malformed equation, so say which step.
    stop("'equation' uses a .patch variable, but 'obj' has no patch table. ",
         "Run a patch-scale primitive first (e.g. msr_area(scale = \"patch\")) ",
         "to decompose the layer into patches.", call. = FALSE)
  } else {
    stop(sprintf(
      "Result length (%d) does not match class (%d), patch (%d), or landscape (1) table.",
      length(result), n_class, n_patch))
  }

  # provenance
  prov <- msk_prov("msr",
                     list(equation = equation, label = label,
                          scale = store_scale, layer = layer))
  obj@provenance <- c(obj@provenance, list(prov))

  obj
}


#' Per-patch row-wise evaluation for matrix-valued patch data
#'
#' For each class, iterates over patches (rows of the matrix). The matrix
#' variable resolves to the focal patch's row (its distances to the other
#' patches); other \code{.patch} variables resolve to the full vector for that
#' class, aligned by patch index. The \code{Inf} diagonal excludes self-terms.
#' The equation must reduce to a scalar per patch.
#'
#' @param env named list of resolved equation variables (metric.scale -> value).
#' @param eq_parsed the parsed equation expression.
#' @param matrix_patch_vars character vector of variable names holding a
#'   matrix-valued patch primitive (e.g. \code{"distance.patch"}).
#' @param scale_map list with class/patch/landscape tables (\code{patch$class}
#'   maps each patch to its class).
#' @return numeric vector of per-patch results, class-blocks in name order.
#' @noRd
.derive_rowwise <- function(env, eq_parsed, matrix_patch_vars, scale_map) {

  # identify classes from the first matrix variable
  first_mat <- env[[matrix_patch_vars[1]]]
  class_names <- names(first_mat)

  # flat patch data uses $class to map patches to classes
  patch_classes <- scale_map$patch$class

  all_results <- numeric(0)

  for (cls in class_names) {

    # indices of this class's patches in the flat patch vectors
    cls_idx <- which(patch_classes == as.integer(cls))
    n_patches_cls <- length(cls_idx)

    for (row_i in seq_len(n_patches_cls)) {

      row_env <- list()

      for (vname in names(env)) {
        val <- env[[vname]]
        parts <- strsplit(vname, ".", fixed = TRUE)[[1]]
        v_scale <- parts[2]

        if (vname %in% matrix_patch_vars) {
          # matrix variable: resolve to this patch's row
          mat <- val[[cls]]
          row_env[[vname]] <- mat[row_i, ]
        } else if (v_scale == "patch") {
          # patch vector: pass the full class subset
          row_env[[vname]] <- val[cls_idx]
        } else {
          # class or landscape variable: pass through as-is
          row_env[[vname]] <- val
        }
      }

      res <- eval(eq_parsed, envir = row_env)
      if (length(res) != 1) {
        stop(sprintf(
          "Equation with matrix patch data must reduce to a scalar per patch, got length %d.",
          length(res)))
      }
      all_results <- c(all_results, res)
    }
  }

  all_results
}


#' Per-patch cell-wise evaluation for cell-level layer data
#'
#' Groups cell values from internal layers (e.g. \code{"_distance"}) by patch ID,
#' then evaluates the equation for each patch. The equation must reduce the
#' cell-value vector to a scalar per patch.
#'
#' @param obj mosaik object (needed to read the patch-ID layer).
#' @param env named list of resolved variables.
#' @param eq_parsed parsed equation.
#' @param cell_vars character vector of variable names with \code{.cell} scale.
#' @param scale_map list with class/patch/landscape tables.
#' @param layer character(1) the input layer name (contains patch IDs).
#' @return numeric vector of per-patch results, aligned with \code{@patches} order.
#' @noRd
.derive_cellwise <- function(obj, env, eq_parsed, cell_vars, scale_map, layer) {

  patch_ids <- msk_pull(obj, layer)
  unique_ids <- sort(unique(patch_ids[!is.na(patch_ids)]))

  # map: patch class from @patches$class (if available)
  patch_classes <- scale_map$patch$class

  all_results <- numeric(0)

  for (pid in unique_ids) {
    cells <- which(patch_ids == pid)

    patch_env <- list()

    for (vname in names(env)) {
      val <- env[[vname]]
      parts <- strsplit(vname, ".", fixed = TRUE)[[1]]
      v_scale <- parts[2]

      if (vname %in% cell_vars) {
        # cell variable: extract values for this patch's cells
        patch_env[[vname]] <- val[cells]
      } else if (v_scale == "patch") {
        # patch vector: pass through (equation should not mix .cell and .patch)
        patch_env[[vname]] <- val
      } else {
        # class or landscape: pass through
        patch_env[[vname]] <- val
      }
    }

    res <- eval(eq_parsed, envir = patch_env)
    if (length(res) != 1) {
      stop(sprintf(
        "Equation with .cell data must reduce to a scalar per patch, got length %d for patch %s.",
        length(res), pid))
    }
    all_results <- c(all_results, res)
  }

  all_results
}
