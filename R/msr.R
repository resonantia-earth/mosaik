#' Measure a metric from an equation
#'
#' Compute a landscape metric from an equation over the classes of a layer,
#' and store the result with them. The classes can be land cover classes or
#' groups, such as the patches numbered by \code{\link{mdf_componentise}}: on a
#' layer of patch numbers, each class is a patch. The values in the equation
#' are primitives, measured by the other \code{msr_*} functions, metrics
#' stored by earlier \code{msr()} calls, or layers, so one metric can build on
#' another.
#'
#' @param obj [`mosaik`]\cr the mosaik with the values already measured.
#' @param equation [`character(1)`][character]\cr the equation that computes
#'   the metric. Each variable in it is named \code{<name>.<focus>}, for
#'   example \code{area.self}, optionally followed by \code{_<layer>} to read
#'   another layer than \code{layer}: \code{area.all_cover}. See Details.
#' @param label [`character(1)`][character]\cr the name the result is stored
#'   under. It must not contain \code{.} or \code{_}, because \code{msr()}
#'   reads these as the start of the focus and of the layer when the label is
#'   used as a variable later. It must not be a layer, a reserved name or a
#'   value already stored with \code{layer}.
#' @param layer [`character(1)`][character]\cr the layer whose classes the
#'   equation is evaluated over, and with which the result is stored.
#'   Defaults to the first layer.
#' @return The mosaik with the result added to the class table of \code{layer}
#'   or to its values for the whole layer; both are read with
#'   \code{\link{msk_table}}.
#' @details
#'   \strong{Focus.} The equation is evaluated with one class of \code{layer}
#'   in focus at a time. The focus of a variable says which classes it reads:
#'   \tabular{ll}{
#'     \code{.self}   \tab the class in focus \cr
#'     \code{.others} \tab every other class of the layer \cr
#'     \code{.all}    \tab all classes of the layer
#'   }
#'   With the areas A 10, B 20 and C 5 and A in focus, \code{area.self} is 10,
#'   \code{area.others} is 20, 5 and \code{area.all} is 10, 20, 5. An equation
#'   with \code{.self} or \code{.others} is evaluated once for each class, and
#'   the results go to the class table. This holds also when it contains
#'   \code{.all}: \code{area.self / sum(area.all)} runs once per class, and
#'   \code{sum(area.all)} is the same total in every run. An equation with only
#'   \code{.all} is evaluated once, and the result is stored for the whole
#'   layer: \code{max(area.all) / sum(area.all)}.
#'
#'   Each run gives the value of one class, or of the whole layer, so it must
#'   come to one number. \code{area.self} is one number,
#'   but \code{area.others} and \code{area.all} hold one number per class, so
#'   a function such as \code{sum()} or \code{max()} has to reduce them. An
#'   equation without that reduction, such as \code{area.all * 2}, has no
#'   single result, and \code{msr()} reports it as an error.
#'
#'   \strong{Names.} A name that is a layer reads cells: \code{canopy.self}
#'   holds the values of the layer "canopy" in the cells of the class in focus,
#'   \code{canopy.all} those in all cells of \code{layer}. Any other name reads
#'   a stored value: a primitive (\code{area}, \code{perimeter},
#'   \code{adjacency}, \code{distance}, \code{dissimilarity}), a column of the
#'   class table, or a metric stored by \code{msr()}. A metric
#'   stored for the whole layer is read with \code{.all}. Some names need no
#'   measurement:
#'   \tabular{ll}{
#'     \code{gid}      \tab the code of each class \cr
#'     \code{complete} \tab \code{TRUE} if a class has no cell on the map
#'       border \cr
#'     \code{x}, \code{y} \tab the coordinates of the cell centres, in map
#'       units
#'   }
#'   These, the primitives and \code{val} and \code{colour} of the class table
#'   are reserved: no layer may take them. A layer that is read by its name,
#'   as in \code{canopy.self}, must not contain \code{.} or \code{_} in that
#'   name, because the \code{.} starts the focus and the \code{_} the layer. Names R knows as constants,
#'   such as \code{pi}, are constants.
#'
#'   \strong{Matrices.} A class by class matrix, such as \code{distance}, is
#'   read along the row of the class in focus: \code{distance.self} is the
#'   entry on the diagonal, the distance of the class to itself (\code{Inf});
#'   \code{distance.others} is the rest of the row, its distances to every
#'   other class, so \code{min(distance.others)} is the distance to the nearest
#'   one; \code{distance.all} is the whole row. With only \code{.all}, the
#'   variable is the whole matrix.
#'
#'   \strong{Other layers.} With \code{_<layer>}, a variable reads another
#'   layer. With \code{.all}, it reads all classes of that layer:
#'   \code{sum(area.all_cover)} is the area of the map. With \code{.self} or
#'   \code{.others}, it reads that layer through the cells of the classes in
#'   focus. After \code{msr_area(layer = "cover")}, on a layer of patches,
#'   \code{area.self_cover} holds, for every cell of the patch in focus, the
#'   area of the cover class that cell belongs to.
#'
#'   \strong{The map border.} A class that touches the map border may extend
#'   beyond it, so its values describe only the part on the map. It still
#'   counts, as a neighbour and in every count. To leave such classes out of
#'   a result, filter with \code{complete}: \code{max(area.all[complete.all])}
#'   is the area of the largest patch that lies wholly on the map.
#'
#'   Each metric is computed once, when \code{msr()} is called; if one of its
#'   primitives is measured again later, the metric is not updated. An
#'   equation can contain any R function, including functions you have written
#'   yourself.
#' @family measure
#' @examples
#' m <- msr_area(landscape, layer = "cover")
#'
#' # the share of the landscape in each class
#' m <- msr(m, equation = "area.self / sum(area.all) * 100", label = "pland",
#'          layer = "cover")
#' msk_table(m, layer = "cover")$pland
#'
#' # a metric built from another metric: Shannon diversity, for the layer
#' m <- msr(m, equation = "-sum(pland.all / 100 * log(pland.all / 100))",
#'          label = "shdi", layer = "cover")
#' msk_table(m, layer = "cover")$shdi
#'
#' # the forest patches: the distance of each to its nearest neighbour, and
#' # the share of the map in the largest patch
#' f <- mdf_filter(m, cover == 47, add = "forest") |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
#'   msr_area(layer = "patch") |>
#'   msr_distance(layer = "patch") |>
#'   msr(equation = "min(distance.others)", label = "enn", layer = "patch") |>
#'   msr(equation = "max(area.all) / sum(area.all_cover) * 100", label = "lpi",
#'       layer = "patch")
#' msk_table(f, layer = "patch")$enn
#' msk_table(f, layer = "patch")$lpi
#'
#' # a layer read through cells: the mean canopy height of each forest patch
#' f <- msr(f, equation = "mean(canopy.self)", label = "height", layer = "patch")
#' msk_table(f, layer = "patch")$height
#' @importFrom checkmate assertClass assertCharacter
#' @export

msr <- function(obj = NULL, equation, label, layer = NULL) {

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = equation, len = 1)
  assertCharacter(x = label, len = 1)
  assertCharacter(x = layer, null.ok = TRUE)

  if (is.null(layer)) layer <- names(obj@layers)[1]
  if (!layer %in% names(obj@layers)) {
    stop("layer '", layer, "' not found in 'obj'.", call. = FALSE)
  }

  # the label is only a name; it must stay reachable, so it may not be a
  # layer or a reserved name, and it must not replace a stored value
  if (!grepl("^[^._]+$", label)) {
    stop(sprintf("'label' must be a name without '.' or '_', not '%s'.", label),
         call. = FALSE)
  }
  if (label %in% c(names(obj@layers), .reserved_names)) {
    stop(sprintf("'label' is '%s', which is a %s; choose another.", label,
                 if (label %in% names(obj@layers)) "layer" else "reserved name"),
         call. = FALSE)
  }
  if (label %in% c(names(obj@categories[[layer]]), names(obj@global[[layer]]))) {
    stop(sprintf("'%s' is already stored with layer '%s'; choose another 'label'.",
                 label, layer), call. = FALSE)
  }

  eq_parsed <- parse(text = equation)
  lv <- msk_pull(obj, layer)
  classes <- .classes_of(obj, layer)
  inside <- which(!is.na(lv))

  # each variable becomes one of three kinds: values per cell, values per
  # class of 'layer', or a value that does not depend on the focus
  notation <- "^([^._]+)\\.(self|others|all)(_(.+))?$"
  vars <- list()
  for (v in all.vars(eq_parsed)) {

    # a value R always knows, such as pi, is a constant
    if (exists(v, envir = baseenv(), inherits = FALSE) &&
        !is.function(get(v, envir = baseenv()))) next

    if (!grepl(notation, v)) {
      stop(sprintf(paste0(
        "Variable '%s' must be named <name>.<focus>, with the focus self, ",
        "others or all, optionally followed by _<layer> (e.g. 'area.self', ",
        "'area.all_cover')."), v), call. = FALSE)
    }
    nm    <- sub(notation, "\\1", v)
    focus <- sub(notation, "\\2", v)
    vl    <- sub(notation, "\\4", v)
    if (!nzchar(vl)) vl <- layer
    if (!vl %in% names(obj@layers)) {
      stop(sprintf("Variable '%s' refers to layer '%s', which is not in 'obj'.",
                   v, vl), call. = FALSE)
    }

    # a name that is a layer reads cells
    if (nm %in% c(names(obj@layers), "x", "y")) {
      if (vl != layer) {
        stop(sprintf(paste0(
          "'%s': '%s' is read through the cells of '%s'; a layer takes no ",
          "_layer suffix."), v, nm, layer), call. = FALSE)
      }
      vars[[v]] <- list(kind = "cell", focus = focus,
                        value = .cell_values(obj, nm))
      next
    }

    # otherwise a stored value: per class of its layer, or one for the layer
    val <- .class_value(obj, vl, nm)
    if (is.null(val)) {
      glob <- obj@global[[vl]][[nm]]
      if (is.null(glob)) {
        stop(sprintf("'%s': no value '%s' is stored with layer '%s'.",
                     v, nm, vl), call. = FALSE)
      }
      if (focus != "all") {
        stop(sprintf(paste0(
          "'%s': '%s' is one value for layer '%s'; read it as '%s.all'."),
          v, nm, vl, nm), call. = FALSE)
      }
      vars[[v]] <- list(kind = "fixed", value = glob)
    } else if (vl == layer) {
      vars[[v]] <- list(kind = "class", focus = focus, value = val)
    } else if (focus == "all") {
      vars[[v]] <- list(kind = "fixed", value = val)
    } else if (is.matrix(val)) {
      stop(sprintf(paste0(
        "'%s': '%s' is a class by class matrix of layer '%s', which can only ",
        "be read as a whole, with .all."), v, nm, vl), call. = FALSE)
    } else {
      # another layer through the cells of the classes in focus
      vals_vl <- msk_pull(obj, vl)
      vars[[v]] <- list(kind = "cell", focus = focus,
                        value = val[match(vals_vl, .classes_of(obj, vl))])
    }
  }

  per_class <- any(vapply(vars, function(x) !is.null(x$focus) &&
                            x$focus != "all", logical(1)))

  if (per_class) {
    members <- lapply(classes, function(k) inside[lv[inside] == k])
    result <- vapply(seq_along(classes), function(i) {
      env <- lapply(vars, .focus_value, i = i, members = members,
                    inside = inside)
      res <- eval(eq_parsed, envir = env)
      if (length(res) != 1) {
        stop(sprintf(paste0(
          "The equation gives %d values for class %s, where one is stored. ",
          "Reduce the values of .others or .all with a function such as ",
          "sum() or max()."), length(res), classes[i]), call. = FALSE)
      }
      as.numeric(res)
    }, numeric(1))
    obj <- .store_class(obj, layer, label, classes, result)
  } else {
    env <- lapply(vars, .focus_value, i = NULL, members = NULL, inside = inside)
    result <- eval(eq_parsed, envir = env)
    if (length(result) != 1) {
      stop(sprintf(paste0(
        "The equation gives %d values for the layer, where one is stored. ",
        "Reduce the values of .all with a function such as sum() or max()."),
        length(result)), call. = FALSE)
    }
    glob <- obj@global
    if (is.null(glob[[layer]])) glob[[layer]] <- list()
    glob[[layer]][[label]] <- result
    obj@global <- glob
  }

  # provenance
  obj <- .update_mosaik(obj, step = step)

  obj
}


# names no layer and no label may take
.reserved_names <- c("area", "perimeter", "adjacency", "distance",
                     "dissimilarity", "gid", "complete", "val", "colour",
                     "x", "y")


#' The classes of a layer
#'
#' From its class table, or else from its values.
#' @noRd
.classes_of <- function(obj, layer) {
  gid <- obj@categories[[layer]]$gid
  if (!is.null(gid)) return(gid)
  v <- msk_pull(obj, layer)
  as.integer(sort(unique(v[!is.na(v)])))
}


#' A value per class of a layer, in the order of its classes
#'
#' A column of the class table, or one of the values that need no
#' measurement: \code{gid} and \code{complete}. \code{NULL} if there is none.
#' @noRd
.class_value <- function(obj, layer, name) {

  classes <- .classes_of(obj, layer)
  if (name == "gid") return(classes)
  if (name == "complete") {
    grid <- matrix(msk_pull(obj, layer), nrow = obj@dims[2],
                   ncol = obj@dims[1], byrow = TRUE)
    border <- c(grid[1, ], grid[nrow(grid), ], grid[, 1], grid[, ncol(grid)])
    return(!classes %in% border)
  }
  obj@categories[[layer]][[name]]
}


#' The values of a layer in every cell
#'
#' A layer, or the coordinates of the cell centres (\code{x}, \code{y}), in
#' map units; cells are numbered row by row from the top-left corner.
#' @noRd
.cell_values <- function(obj, name) {

  if (!name %in% c("x", "y")) return(msk_pull(obj, name))
  ext <- obj@extent
  res <- msk_res(obj)
  d <- obj@dims
  if (name == "x") {
    rep(ext[1] + (seq_len(d[1]) - 0.5) * res[1], times = d[2])
  } else {
    rep(ext[4] - (seq_len(d[2]) - 0.5) * res[2], each = d[1])
  }
}


#' A variable's value with class i in focus
#'
#' @param var list: \code{kind} (cell, class or fixed), \code{focus} and
#'   \code{value}.
#' @param i the position of the class in focus, \code{NULL} when the equation
#'   is evaluated once for the layer.
#' @param members the cells of each class.
#' @param inside the cells of the layer that are not \code{NA}.
#' @noRd
.focus_value <- function(var, i, members, inside) {

  val <- var$value
  switch(var$kind,
    fixed = val,
    cell = if (is.null(i)) val[inside] else switch(var$focus,
      self   = val[members[[i]]],
      others = val[setdiff(inside, members[[i]])],
      all    = val[inside]),
    class = if (is.null(i)) val else if (is.matrix(val)) switch(var$focus,
      self   = val[i, i],
      others = val[i, -i],
      all    = val[i, ]) else switch(var$focus,
      self   = val[i],
      others = val[-i],
      all    = val))
}
