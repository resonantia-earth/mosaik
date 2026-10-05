#' Measure a composite metric
#'
#' Compute a metric from an equation over values that \code{msr_*} functions
#' or earlier \code{msr()} calls have stored. Variables in the equation use
#' \code{metric.scale} notation to say which table each value comes from.
#'
#' @param obj [`mosaik`]\cr the mosaik with the values already measured.
#' @param equation [`character(1)`][character]\cr a mathematical expression using
#'   \code{metric.scale} notation, where scale is one of \code{class},
#'   \code{patch}, \code{landscape}, or \code{cell}, optionally followed by
#'   \code{_layer}. See Details.
#' @param label [`character(1)`][character]\cr the name of the metric, written
#'   like a variable: \code{"pland.class"} stores the result with the classes,
#'   \code{"enn.patch"} with the patches, \code{"shdi.landscape"} as a
#'   landscape value. The name before the \code{.} must not exist at that level
#'   already. Without the level (\code{"pland"}), the result's length decides
#'   where it goes; see Details.
#' @param layer [`character(1)`][character]\cr the layer whose results the
#'   metric is stored with, and whose values a variable without \code{_layer}
#'   refers to. Defaults to the first layer.
#' @return The mosaik with the metric added to the results of \code{layer}.
#' @details
#'   \strong{Variables.} A variable can name any value stored at its level: a
#'   primitive, or a metric computed by an earlier \code{msr()} call, so a
#'   metric can be built from other metrics. Such a metric is computed once,
#'   when \code{msr()} runs; measuring a primitive again afterwards does not
#'   update it. Names R always knows, such as \code{pi}, are constants, and any
#'   function R can find, including your own, can be used.
#'
#'   \strong{Where the result is stored.} With the level in the label,
#'   \code{msr()} stores the result there and stops if the result does not
#'   fit: one value per class, one per patch, or a single value for the
#'   landscape. Without it, a result with one value per class goes to the
#'   classes, one per patch to the patches, a single value to the landscape.
#'   That guess fails where the counts coincide: on a layer with one class, a
#'   landscape value such as Shannon diversity would be stored as a class
#'   value. Name the level in the label to be sure.
#'
#'   \strong{Several layers.} \code{area.class} refers to \code{layer};
#'   \code{area.class_core} to the layer \code{core}. The metric name contains
#'   no \code{.} and the scale no \code{_}, so everything after the first
#'   \code{_} following the scale is the layer name, underscores included.
#'   Landscape values of different layers combine freely. Class values of
#'   another layer combine with those of \code{layer} if both layers have the
#'   same classes, as two binary layers do; to use a total over the classes of
#'   another layer, compute it on that layer first (a landscape value) and
#'   refer to that. Patch values of another layer only combine if both layers
#'   have the same patches; relating patches of different layers is a job for
#'   \code{\link{mdf_summarise}}.
#'
#'   \strong{Patches and cells.} Patches are those numbered by
#'   \code{\link{mdf_componentise}} on \code{layer}; a \code{.patch} or
#'   \code{.cell} variable needs that step first. The \code{cell} scale reads
#'   the per-cell layer \code{"_<metric>_<layer>"} (e.g.
#'   \code{distance.cell_forest} reads \code{"_distance_forest"}, written by
#'   \code{msr_cost(scale = "cell", layer = "forest")}). The equation is then
#'   evaluated once per patch, on the values of that patch's cells, and must
#'   give one number per patch.
#'
#'   \strong{Distances between patches.} \code{distance.patch} holds, for every
#'   class, a matrix of the distances between its patches. An equation that
#'   uses it is evaluated once per patch: \code{distance.patch} is that patch's
#'   row, its distances to the other patches of its class. Any other
#'   \code{.patch} variable holds the values of those same patches, in the same
#'   order, so \code{sum(area.patch / distance.patch^2)} adds up the area of
#'   every other patch divided by its squared distance (the proximity index).
#'   A patch's distance to itself is \code{Inf}, so it drops out of such sums
#'   on its own. The equation must give one number per patch; so
#'   \code{"min(distance.patch)"} is each patch's nearest-neighbour distance.
#'
#'   Costs are named after what they measure: \code{distance} without a cost
#'   surface, otherwise the name of the cost layer, e.g. \code{friction.patch}
#'   (see \code{\link{msr_cost}}).
#' @family measure
#' @examples
#' m <- landscape |>
#'   msr_area(scale = "class", layer = "cover") |>
#'   msr_area(scale = "landscape", layer = "cover")
#'
#' # the share of the landscape in each class
#' m <- msr(m, equation = "area.class / area.landscape * 100",
#'          label = "pland.class", layer = "cover")
#' msk_categories(m, layer = "cover")$pland
#'
#' # a metric built from another metric: Shannon diversity
#' m <- msr(m, equation = "-sum(pland.class / 100 * log(pland.class / 100))",
#'          label = "shdi.landscape", layer = "cover")
#' msk_global(m, layer = "cover")$shdi
#'
#' # combining two layers: the share of the forest that is core
#' f <- mdf_binarise(landscape, match = 47, layer = "cover", add = "forest") |>
#'   mdf_erode(layer = "forest", add = "core") |>
#'   msr_area(scale = "class", layer = "forest") |>
#'   msr_area(scale = "class", layer = "core") |>
#'   msr(equation = "area.class_core / area.class_forest",
#'       label = "core_share.class", layer = "core")
#' msk_categories(f, layer = "core")$core_share
#'
#' # per patch: the distance of each forest patch to its nearest neighbour
#' f <- f |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
#'   msr_cost(scale = "patch", layer = "forest") |>
#'   msr(equation = "min(distance.patch)", label = "enn.patch",
#'       layer = "forest")
#' msk_patches(f, layer = "forest")$enn
#' @importFrom checkmate assertCharacter assertClass
#' @export

msr <- function(obj = NULL, equation, label, layer = NULL) {

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = equation, len = 1)
  assertCharacter(x = label, len = 1)
  assertCharacter(x = layer, null.ok = TRUE)

  if (is.null(layer)) layer <- names(obj@layers)[1]

  # the label is a name, optionally followed by the level the result is
  # stored at, in the notation of the variables
  label_notation <- "^([^.]+)(\\.(class|patch|landscape))?$"
  if (!grepl(label_notation, label)) {
    stop(sprintf(paste0(
      "'label' must be a name without '.', optionally followed by .class, ",
      ".patch or .landscape (e.g. 'pland.class'), not '%s'."), label),
      call. = FALSE)
  }
  label_scale <- sub(label_notation, "\\3", label)
  label <- sub(label_notation, "\\1", label)

  # the results of one layer at each level
  tables_of <- function(l) {
    list(class     = obj@categories[[l]],
         patch     = obj@patches[[l]],
         landscape = obj@global[[l]])
  }
  scale_map <- tables_of(layer)

  # parse equation variables
  eq_parsed <- parse(text = equation)
  eq_vars <- all.vars(eq_parsed)

  # resolve each variable via metric.scale notation, optionally followed by
  # _layer: the metric has no '.', the scale no '_', so the split is unique
  env <- list()
  var_scale <- character(0)
  matrix_patch_vars <- character(0)
  cell_vars <- character(0)
  notation <- "^([^.]+)\\.(class|patch|landscape|cell)(_(.+))?$"

  for (v in eq_vars) {

    # a value R always knows, such as pi, is a constant
    if (exists(v, envir = baseenv(), inherits = FALSE) &&
        !is.function(get(v, envir = baseenv()))) next

    if (!grepl(notation, v)) {
      stop(sprintf(paste0(
        "Variable '%s' must use metric.scale notation, optionally followed by ",
        "_layer (e.g. 'area.class', 'area.class_forest', 'distance.cell')."), v),
        call. = FALSE)
    }

    metric  <- sub(notation, "\\1", v)
    v_scale <- sub(notation, "\\2", v)
    v_layer <- sub(notation, "\\4", v)
    if (!nzchar(v_layer)) v_layer <- layer
    if (!v_layer %in% names(obj@layers)) {
      stop(sprintf("Variable '%s' refers to layer '%s', which is not in 'obj'.",
                   v, v_layer), call. = FALSE)
    }
    var_scale[v] <- v_scale

    if (v_scale == "cell") {
      # cell scale: look up the internal layer "_<metric>_<layer>"
      cell_layer <- paste0("_", metric, "_", v_layer)
      if (!(cell_layer %in% names(obj@layers))) {
        stop(sprintf("Internal layer '%s' not found. Run a measure at scale = 'cell' on layer '%s' first.",
                     cell_layer, v_layer), call. = FALSE)
      }
      env[[v]] <- msk_pull(obj, cell_layer)
      cell_vars <- c(cell_vars, v)
      next
    }

    # patches of two layers only correspond if both were numbered into the
    # same patches, in the same order
    same_patches <- function(a, b) {
      pa <- obj@patches[[a]]
      pb <- obj@patches[[b]]
      !is.null(pa$ids) && !is.null(pb$ids) &&
        all(c(pa$ids, pb$ids) %in% names(obj@layers)) &&
        identical(msk_pull(obj, pa$ids), msk_pull(obj, pb$ids)) &&
        identical(pa$patch, pb$patch)
    }
    if (v_scale == "patch" && v_layer != layer && !same_patches(v_layer, layer)) {
      stop(sprintf(paste0(
        "'%s' is a patch value of layer '%s', whose patches differ from those ",
        "of '%s'. To relate patches of different layers, summarise one within ",
        "the other with mdf_summarise()."), v, v_layer, layer), call. = FALSE)
    }

    tbl <- tables_of(v_layer)[[v_scale]]
    if (is.null(tbl) || is.null(tbl[[metric]])) {
      stop(sprintf("Metric '%s' not found at %s level for layer '%s'.",
                   metric, v_scale, v_layer), call. = FALSE)
    }

    # class values are combined position by position, so both layers must
    # list the same classes in the same order
    if (v_scale == "class" && v_layer != layer &&
        !identical(as.numeric(tbl$gid), as.numeric(scale_map$class$gid))) {
      stop(sprintf(paste0(
        "'%s' is a class value of layer '%s', whose classes (%s) differ from ",
        "those of '%s' (%s)."), v, v_layer, paste(tbl$gid, collapse = ", "),
        layer, paste(scale_map$class$gid, collapse = ", ")), call. = FALSE)
    }
    val <- tbl[[metric]]

    # detect matrix-valued patch data (e.g. distance stored as list of matrices)
    if (v_scale == "patch" && is.list(val) && !is.data.frame(val)) {
      matrix_patch_vars <- c(matrix_patch_vars, v)
    }

    env[[v]] <- val
  }

  has_cell <- length(cell_vars) > 0
  if (has_cell) {
    result <- .derive_cellwise(env, eq_parsed, cell_vars,
                               .patches_of(obj, layer))
  } else if (length(matrix_patch_vars) > 0) {
    result <- .derive_rowwise(env, eq_parsed, matrix_patch_vars, var_scale,
                              scale_map$patch$class)
  } else {
    result <- eval(eq_parsed, envir = env)
  }

  # the level the result is stored at: as the label says, which the result
  # must fit; without one, from the result's length
  n_class <- length(scale_map$class$gid)
  n_patch <- length(scale_map$patch$patch)

  if (nzchar(label_scale)) {
    if (label_scale == "patch" && n_patch == 0) {
      stop("'label' stores the result with the patches, but layer '", layer,
           "' has no patches. Number them first with mdf_componentise(layer = \"",
           layer, "\", ...).", call. = FALSE)
    }
    n_expected <- c(class = n_class, patch = n_patch, landscape = 1L)[[label_scale]]
    if (length(result) != n_expected) {
      takes <- switch(label_scale,
        class = sprintf("one value per class of layer '%s' (%d)", layer, n_class),
        patch = sprintf("one value per patch of layer '%s' (%d)", layer, n_patch),
        landscape = "one value")
      stop(sprintf("'%s.%s' takes %s, but the equation gives %d.",
                   label, label_scale, takes, length(result)), call. = FALSE)
    }
    store_scale <- label_scale
  } else if (has_cell) {
    # cell-wise evaluation always produces per-patch results
    store_scale <- "patch"
  } else if (length(result) == n_class && n_class > 0) {
    store_scale <- "class"
  } else if (length(result) == n_patch && n_patch > 0) {
    store_scale <- "patch"
  } else if (length(result) == 1) {
    store_scale <- "landscape"
  } else if (n_patch == 0 && any(var_scale == "patch")) {
    # the equation asks for patch scale but the layer has no patches yet, so
    # there is no table to store the result in. This is a missing step rather
    # than a malformed equation, so say which step.
    stop("'equation' uses a .patch variable, but layer '", layer, "' has no ",
         "patches. Number them first with mdf_componentise(layer = \"", layer,
         "\", ...).", call. = FALSE)
  } else {
    stop(sprintf(
      "Result length (%d) does not match class (%d), patch (%d), or landscape (1) table.",
      length(result), n_class, n_patch))
  }

  # a label must not replace a stored value: that would silently overwrite a
  # primitive, an earlier metric or the identifiers the table is keyed by
  if (label %in% names(scale_map[[store_scale]])) {
    stop(sprintf("'%s' already exists at %s level for layer '%s'; choose another 'label'.",
                 label, store_scale, layer), call. = FALSE)
  }
  slot_name <- c(class = "categories", patch = "patches",
                 landscape = "global")[[store_scale]]
  if (is.null(methods::slot(obj, slot_name)[[layer]])) {
    tbl <- methods::slot(obj, slot_name)
    tbl[[layer]] <- list()
    methods::slot(obj, slot_name) <- tbl
  }
  switch(store_scale,
         class     = obj@categories[[layer]][[label]] <- result,
         patch     = obj@patches[[layer]][[label]] <- result,
         landscape = obj@global[[layer]][[label]] <- result)

  # provenance
  obj <- .update_mosaik(obj, step = step)

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
#' @param var_scale named character, the scale of each variable.
#' @param patch_classes the class of each patch of the layer, in roster order.
#' @return numeric vector of per-patch results, in the order of the patch
#'   record.
#' @noRd
.derive_rowwise <- function(env, eq_parsed, matrix_patch_vars, var_scale,
                            patch_classes) {

  # identify classes from the first matrix variable
  first_mat <- env[[matrix_patch_vars[1]]]
  class_names <- names(first_mat)

  # one result per patch, at that patch's position in the record
  all_results <- rep(NA_real_, length(patch_classes))

  for (cls in class_names) {

    # indices of this class's patches in the flat patch vectors
    cls_idx <- which(patch_classes == as.integer(cls))
    n_patches_cls <- length(cls_idx)

    for (row_i in seq_len(n_patches_cls)) {

      row_env <- list()

      for (vname in names(env)) {
        val <- env[[vname]]
        v_scale <- var_scale[[vname]]

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
      all_results[cls_idx[row_i]] <- res
    }
  }

  all_results
}


#' Per-patch cell-wise evaluation for cell-level layer data
#'
#' Groups cell values from internal layers (e.g. \code{"_distance_forest"}) by
#' patch, then evaluates the equation for each patch. The equation must reduce
#' the cell-value vector to a scalar per patch.
#'
#' @param env named list of resolved variables.
#' @param eq_parsed parsed equation.
#' @param cell_vars character vector of variable names with \code{.cell} scale.
#' @param patches the patches of the layer, from \code{.patches_of()}.
#' @return numeric vector of per-patch results, in the order of the patch
#'   record.
#' @noRd
.derive_cellwise <- function(env, eq_parsed, cell_vars, patches) {

  patch_ids <- patches$ids

  all_results <- numeric(0)

  for (pid in patches$patch) {
    cells <- which(!is.na(patch_ids) & patch_ids == pid)

    patch_env <- list()

    for (vname in names(env)) {
      val <- env[[vname]]
      # a cell variable becomes this patch's cells; anything else (class,
      # landscape, patch vector) passes through unchanged
      patch_env[[vname]] <- if (vname %in% cell_vars) val[cells] else val
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
