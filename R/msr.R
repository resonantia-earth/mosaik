#' Measure a metric from an equation
#'
#' Compute a landscape metric from an equation over values that have already
#' been measured, and store the result with them. The values are primitives,
#' measured by the other \code{msr_*} functions, or metrics stored by earlier
#' \code{msr()} calls, so one metric can build on another.
#'
#' @param obj [`mosaik`]\cr the mosaik with the values already measured.
#' @param equation [`character(1)`][character]\cr the equation that computes
#'   the metric. Each variable in it refers to a measured value and is named
#'   \code{<metric>.<scale>}, for example \code{area.class} (the area of each
#'   class) or \code{perimeter.patch} (the perimeter of each patch). The scale
#'   is one of class, patch or landscape. A variable refers to the values of
#'   \code{layer}; to use the values of another layer, add that layer's name
#'   after an underscore: \code{area.landscape_cover} is the area of the layer
#'   "cover". A variable with the scale cell names a layer instead:
#'   \code{canopy.cell} holds the values of the layer "canopy" in every cell.
#' @param label [`character(1)`][character]\cr the name the result is stored
#'   under, in the same notation as the variables, for example
#'   \code{"pland.class"}, \code{"enn.patch"} or \code{"shdi.landscape"}. The
#'   scale of the label names the group the result belongs to: each class of
#'   \code{layer}, each patch, or the whole layer. The equation is evaluated
#'   once for each group, and every variable in it holds the values of that
#'   group. For example, \code{canopy.cell} then holds the canopy height of the
#'   cells in that group. The equation must return one value for each group;
#'   otherwise \code{msr()} stops with an error. A label without a scale
#'   (\code{"pland"}) is stored according to the number of values the equation
#'   returns. The label must not already exist at that scale of \code{layer}.
#' @param layer [`character(1)`][character]\cr the layer whose values the
#'   variables refer to by default, and with which the result is stored.
#'   Defaults to the first layer.
#' @return The mosaik with the result added to the results of \code{layer}.
#' @details
#'   \strong{Variables.} An equation can use the primitives, which the other
#'   \code{msr_*} functions measure (area, perimeter, number, adjacency, cost
#'   and dissimilarity), and any metric that an earlier \code{msr()} call has
#'   stored. In this way, a metric can be built from other metrics. Each
#'   metric is computed once, when \code{msr()} is called; if one of its
#'   primitives is measured again later, the metric is not updated. An
#'   equation can also contain constants such as \code{pi} and any R function,
#'   including functions you have written yourself.
#'
#'   \strong{Combining layers.} Landscape values of different layers can
#'   always be combined, because each layer has exactly one. Class values of
#'   two layers can only be combined if both layers have the same classes, as
#'   two binary layers do. Patch values of two layers can only be combined if
#'   both layers were numbered into the same patches. If they were not, for
#'   example the forest patches and the core areas within them, use the other
#'   layer as a cell variable: \code{mean(core.cell)} with the label
#'   \code{"coreshare.patch"} gives the share of core in each forest patch, if
#'   the layer "core" is 1 for core and 0 for every other cell.
#'
#'   \strong{Patches.} Patch values require that the patches of \code{layer}
#'   have been numbered with \code{\link{mdf_componentise}}. Patch metrics often
#'   concern only some of the classes, for example those that make up a
#'   habitat. In that case, combine these classes into one binary layer with
#'   \code{\link{mdf_binarise}}, number the patches of that layer and compute
#'   the metric on it. A patch that touches the map border has \code{NA} for
#'   all its values, because the map shows only part of it. To leave such
#'   patches out of a sum or a maximum, use \code{na.rm = TRUE}.
#'
#'   \strong{Cell values.} Every layer can be used as a cell variable: land
#'   cover, canopy height, or a layer made by an \code{mdf_*} function, such as
#'   the distance of every cell to the edge of its patch from
#'   \code{\link{mdf_distance}}. The layer name must be a valid R name, so it
#'   must not contain a dot or start with an underscore. A label without a
#'   scale groups an equation with cell variables by patch.
#'
#'   \strong{Distances between patches.} \code{\link{msr_distance}} stores the
#'   distances between the patches of each class as a matrix. Within a patch,
#'   \code{distance.patch} holds the distances from this patch to all other
#'   patches of the same class. Thus \code{min(distance.patch)} is the
#'   distance to the nearest patch. The distance of a patch to itself is
#'   \code{Inf}, so it does not count. A label without a scale groups an
#'   equation with \code{distance.patch} by patch. The matrix is called
#'   \code{distance}, unless \code{msr_distance()} was given a layer with the
#'   cost of crossing each cell (its argument \code{cost}). Then the matrix
#'   takes the name of that layer: after \code{msr_distance(cost = "friction")},
#'   the variable is \code{friction.patch}.
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
#' # combining two layers: the share of the forest that is core, from the
#' # forest class (1) of both layers
#' f <- mdf_filter(landscape, cover == 47, add = "forest") |>
#'   mdf_erode(layer = "forest", add = "core") |>
#'   msr_area(scale = "class", layer = "forest") |>
#'   msr_area(scale = "class", layer = "core") |>
#'   msr(equation = "area.class_core[gid.class == 1] / area.class[gid.class == 1]",
#'       label = "core_share.landscape", layer = "forest")
#' msk_global(f, layer = "forest")$core_share
#'
#' # per patch: the distance of each forest patch to its nearest neighbour;
#' # patches the map border cuts are NA
#' f <- f |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
#'   msr_distance(layer = "forest") |>
#'   msr(equation = "min(distance.patch)", label = "enn.patch",
#'       layer = "forest")
#' msk_patches(f, layer = "forest")$enn
#'
#' # a layer as cell variable: the mean canopy height of each forest patch,
#' # and of each land cover class
#' f <- msr(f, equation = "mean(canopy.cell)", label = "height.patch",
#'          layer = "forest")
#' msk_patches(f, layer = "forest")$height
#' f <- msr(f, equation = "mean(canopy.cell)", label = "height.class",
#'          layer = "cover")
#' msk_categories(f, layer = "cover")$height
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
        "_layer (e.g. 'area.class', 'area.class_forest'), or name a layer as ",
        "<layer>.cell (e.g. 'canopy.cell')."), v),
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
      # a cell variable names a layer: its values in every cell
      if (sub(notation, "\\4", v) != "") {
        stop(sprintf(paste0(
          "'%s': a cell variable names a layer itself and takes no _layer ",
          "suffix."), v), call. = FALSE)
      }
      if (!(metric %in% names(obj@layers))) {
        stop(sprintf("Variable '%s' refers to layer '%s', which is not in 'obj'.",
                     v, metric), call. = FALSE)
      }
      env[[v]] <- msk_pull(obj, metric)
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
  n_patch <- length(scale_map$patch$patch)

  # the group the result belongs to: as the label says; a label without one
  # groups by patch if the equation reads cells or distances, otherwise not
  group <- if (nzchar(label_scale)) {
    label_scale
  } else if (has_cell || length(matrix_patch_vars) > 0) {
    "patch"
  } else {
    ""
  }

  if (group == "patch" && n_patch == 0) {
    stop("'label' stores the result with the patches, but layer '", layer,
         "' has no patches. Number them first with mdf_componentise(layer = \"",
         layer, "\", ...).", call. = FALSE)
  }

  # the classes of the layer, from its class table or else from its values
  classes <- scale_map$class$gid
  if (group == "class" && is.null(classes)) {
    v <- msk_pull(obj, layer)
    classes <- sort(unique(v[!is.na(v)]))
  }

  if (nzchar(group)) {
    result <- .derive_grouped(obj, layer, env, eq_parsed, var_scale, cell_vars,
                              matrix_patch_vars, group, classes)
    store_scale <- group
  } else {
    result <- eval(eq_parsed, envir = env)
  }

  # without a level in the label, the result's length decides
  if (!nzchar(group)) {
    n_class <- length(scale_map$class$gid)
    if (length(result) == n_class && n_class > 0) {
      store_scale <- "class"
    } else if (length(result) == n_patch && n_patch > 0) {
      store_scale <- "patch"
    } else if (length(result) == 1) {
      store_scale <- "landscape"
    } else if (n_patch == 0 && any(var_scale == "patch")) {
      # the equation asks for patch scale but the layer has no patches yet, so
      # there is no table to store the result in. This is a missing step
      # rather than a malformed equation, so say which step.
      stop("'equation' uses a .patch variable, but layer '", layer, "' has no ",
           "patches. Number them first with mdf_componentise(layer = \"", layer,
           "\", ...).", call. = FALSE)
    } else {
      stop(sprintf(
        "Result length (%d) does not match class (%d), patch (%d), or landscape (1) table.",
        length(result), n_class, n_patch))
    }
  }

  # a class result on a layer without class table starts one
  if (store_scale == "class" && is.null(scale_map$class$gid)) {
    obj@categories[[layer]]$gid <- classes
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


#' Evaluate an equation once for every group
#'
#' The group is the scale of the label: each patch, each class of the layer,
#' or the whole layer. In every group, each variable holds the values of that
#' group: a cell variable the values of the group's cells, a patch or class
#' variable the value of the group (for a patch, also the value of its class;
#' for a class, the values of its patches), a distance matrix the patch's own
#' row. Landscape values are the same in every group.
#'
#' @param obj the mosaik.
#' @param layer the layer the result is stored with.
#' @param env named list of resolved variables.
#' @param eq_parsed parsed equation.
#' @param var_scale named character, the scale of each variable.
#' @param cell_vars names of the cell variables.
#' @param matrix_vars names of the variables holding distance matrices.
#' @param group \code{"patch"}, \code{"class"} or \code{"landscape"}.
#' @param classes the classes of \code{layer}.
#' @return one value per group, in the order of the patch record or the
#'   classes.
#' @noRd
.derive_grouped <- function(obj, layer, env, eq_parsed, var_scale, cell_vars,
                            matrix_vars, group, classes) {

  if (group == "landscape") {
    # every cell of the layer, every value as it is
    if (length(cell_vars)) {
      cells <- which(!is.na(msk_pull(obj, layer)))
      for (v in cell_vars) env[[v]] <- env[[v]][cells]
    }
    res <- eval(eq_parsed, envir = env)
    if (length(res) != 1) {
      stop(sprintf("'.landscape' takes one value, but the equation gives %d.",
                   length(res)), call. = FALSE)
    }
    return(res)
  }

  rec <- obj@patches[[layer]]
  patch_class <- rec$class
  linked <- !is.null(rec$ids) && rec$ids %in% names(obj@layers)
  patches <- if (linked) .patches_of(obj, layer) else NULL

  # the cells of every group
  if (length(cell_vars)) {
    if (group == "patch") {
      if (!linked) .patches_of(obj, layer)   # stops with the reason
      ids <- patches$ids
      members <- lapply(rec$patch, function(p) which(!is.na(ids) & ids == p))
    } else {
      vals <- msk_pull(obj, layer)
      members <- lapply(classes, function(k) which(!is.na(vals) & vals == k))
    }
  }

  n <- if (group == "patch") length(rec$patch) else length(classes)
  result <- rep(NA_real_, n)

  for (i in seq_len(n)) {

    k <- if (group == "patch") patch_class[i] else classes[i]
    g_env <- list()

    for (v in names(env)) {
      val <- env[[v]]
      g_env[[v]] <- if (v %in% cell_vars) {
        val[members[[i]]]
      } else if (v %in% matrix_vars) {
        # a patch's own distances: its row in the matrix of its class
        mat <- val[[as.character(k)]]
        if (group == "patch") mat[as.character(rec$patch[i]), ] else mat
      } else if (var_scale[[v]] == "patch") {
        if (group == "patch") val[i] else val[patch_class == k]
      } else if (var_scale[[v]] == "class") {
        val[match(k, classes)]
      } else {
        val
      }
    }

    res <- eval(eq_parsed, envir = g_env)
    if (length(res) != 1) {
      stop(sprintf(paste0(
        "The equation must give one value for each %s, but gives %d for %s %s."),
        group, length(res), group, if (group == "patch") rec$patch[i] else k),
        call. = FALSE)
    }
    result[i] <- res
  }

  # a patch the map border cuts has NA for all its values
  if (group == "patch" && linked) result[patches$clipped] <- NA

  result
}
