#' Capture the call of a mosaik function
#'
#' Called on the first line of every function that changes a mosaik, before
#' the function touches its own arguments. It records what the function was
#' asked to do, so that \code{.update_mosaik()} can later write the provenance
#' entry without the function stating it.
#'
#' @details
#' \strong{Arguments.} Every formal argument except \code{obj} is recorded with
#' its value, defaults included. Values are read from the calling frame, which
#' forces the argument's promise once, so an argument such as
#' \code{runif(10)} is not evaluated a second time. Arguments listed in
#' \code{.quoted_args} are expressions over layers of the mosaik and are kept
#' unevaluated. Values larger than \code{.max_arg_size} are stored as a digest
#' under \code{used} rather than copied into the history.
#'
#' \strong{Nesting.} When a mosaik function is called by another one that
#' records itself (\code{mdf_interpolate} calling \code{mdf_morph},
#' \code{mdf_loop} replaying its recipe), only the outer call is recorded.
#' \code{mdf()} records nothing itself, so each step it replays is recorded.
#'
#' \strong{Layers.} The input layer (\code{layer}, or \code{layers}) and the
#' output layer (\code{add}, else the input; \code{name} for \code{drw_*}) are
#' resolved here, so that the entry names them and the function need not.
#'
#' @return an object of class \code{msk_step}.
#' @noRd

.step <- function() {

  frame <- parent.frame()
  fun   <- sys.function(-1)
  call  <- sys.call(-1)
  outer <- parent.frame(2)
  ns    <- topenv(environment(fun))
  fn    <- .fn_name(call[[1]], fun, ns)

  # an enclosing mosaik function already records this call
  above <- utils::head(sys.frames(), -2)
  nested <- any(vapply(above, function(f) {
    inherits(get0("step", envir = f, inherits = FALSE), "msk_step")
  }, logical(1)))

  mc <- match.call(fun, call, expand.dots = FALSE, envir = outer)
  quoted <- .quoted_args[[fn]]
  is_missing <- function(a) eval(call("missing", as.name(a)), frame)

  args <- list()
  used <- list()
  keep_value <- function(nm, v) {
    if (utils::object.size(v) > .max_arg_size && !.is_recipe_shell(v)) {
      d <- digest::digest(v)
      used[[nm]] <<- list(class = class(v)[1], digest = d)
      v <- structure(list(digest = d), class = "msk_used")
    }
    args[nm] <<- list(v)
  }

  for (a in setdiff(names(formals(fun)), "obj")) {
    if (a == "...") {
      dots <- mc$...
      nms <- names(dots)
      if (is.null(nms)) nms <- rep("", length(dots))
      for (i in seq_along(dots)) {
        # an unnamed symbol names a layer; anything else is a value
        if (!nzchar(nms[i]) && is.name(dots[[i]])) {
          args[[length(args) + 1]] <- as.character(dots[[i]])
        } else {
          v <- eval(call("...elt", i), frame)
          if (nzchar(nms[i])) keep_value(nms[i], v) else args[length(args) + 1] <- list(v)
        }
      }
    } else if (a %in% quoted) {
      e <- mc[[a]]
      # when replayed, the expression arrives bound to a symbol; unwrap it
      if (is.name(e) && exists(as.character(e), envir = outer, inherits = FALSE)) {
        bound <- get(as.character(e), envir = outer, inherits = FALSE)
        if (is.call(bound) || is.name(bound)) e <- bound
      }
      args[a] <- list(e)
    } else if (!is_missing(a) || !identical(formals(fun)[[a]], quote(expr = ))) {
      keep_value(a, get(a, envir = frame))
    }
  }

  # input and output layer
  obj <- if ("obj" %in% names(formals(fun)) && !is_missing("obj")) get("obj", envir = frame)
  has_layers <- methods::is(obj, "mosaik") && length(obj@layers) > 0
  from <- if ("layers" %in% names(args)) args$layers else args$layer
  if (!is.character(from)) from <- NULL
  reads_layer <- any(c("layer", "layers") %in% names(formals(fun)))
  if (is.null(from) && has_layers && reads_layer) {
    from <- if ("layers" %in% names(args)) names(obj@layers) else names(obj@layers)[1]
  }
  # a measure's 'name' names the value it stores, not a layer
  to <- if ("name" %in% names(args) && !startsWith(fn, "msr")) {
    args$name
  } else if ("add" %in% names(formals(fun))) {
    if (is.null(args$add)) from[1] else args$add
  }

  # every other layer the step reads: the layers named by other arguments and
  # those a predicate refers to. Without layers (a recipe) they cannot be
  # checked, so every name is kept except the keywords of mdf_distance
  other <- character()
  for (a in intersect(c("by", "source", "origin", "anchor", "cost"), names(args))) {
    if (is.character(args[[a]])) other <- c(other, args[[a]])
  }
  for (a in quoted) other <- c(other, all.vars(args[[a]]))
  other <- if (has_layers) intersect(other, names(obj@layers)) else
    setdiff(other, c("foreground", "background"))
  # a mask written by mdf_filter does not read the layer it is written into
  if (fn == "mdf_filter" && !isTRUE(args$value)) from <- NULL
  from <- unique(c(from, other))
  if (length(from) == 0) from <- NULL

  structure(list(fn = fn, args = args, used = used, from = from, to = to,
                 agent = .agent(ns), record = !nested),
            class = "msk_step")
}


#' Write to a mosaik and record the step
#'
#' The one place where a mosaik is changed. Replaces the slots given in
#' \code{...}, writes \code{values} as the output layer of \code{step},
#' compresses new layers, appends the provenance entry and validates the
#' result.
#'
#' @param obj the mosaik to change; \code{NULL} starts a recipe.
#' @param ... slots to replace, by name (\code{layers}, \code{categories},
#'   \code{global}, \code{extent}, \code{dims}, \code{crs}).
#' @param values cell values of the output layer \code{step$to}.
#' @param gid,val for a categorical output layer: its group IDs and labels.
#'   Only those that occur in \code{values} are registered.
#' @param keep when overwriting a continuous layer in place, keep the fields
#'   other packages attached to its category entry. \code{FALSE} for an
#'   operator that changes what the layer is (binarise, componentise).
#' @param step the \code{msk_step} captured on the function's first line, or
#'   \code{NULL} for a write that is not a user step (a cache layer).
#' @param activity further PROV verbs for \code{wasGeneratedBy}, for what the
#'   function found out while running (\code{mdf_loop}'s iterations).
#' @return the changed mosaik.
#' @noRd

.update_mosaik <- function(obj, ..., values = NULL, gid = NULL, val = NULL,
                           keep = TRUE, step = NULL, activity = list()) {

  if (is.null(obj)) obj <- .recipe_shell()

  slots <- list(...)
  if (!is.null(slots$dims)) slots$dims <- as.integer(slots$dims)
  for (s in names(slots)) methods::slot(obj, s) <- slots[[s]]
  changed <- if ("layers" %in% names(slots)) names(obj@layers) else character()

  if (!is.null(values)) {
    layer <- step$to
    prior <- obj@categories[[layer]]
    # a class table describes the old values, so only the fields of a layer
    # without classes (attached by another package) survive the rewrite
    extra <- if (keep && is.null(gid) && !is.null(prior) && is.null(prior$gid)) {
      prior[setdiff(names(prior), c("gid", "val"))]
    } else list()

    obj@layers[[layer]] <- values
    obj@categories[[layer]] <- NULL
    # what was measured on the old values no longer describes the layer
    obj@global[[layer]] <- NULL
    if (!is.null(gid)) {
      present <- gid %in% unique(values[!is.na(values)])
      obj@categories[[layer]] <- list(gid = gid[present], val = val[present])
    } else if (length(extra)) {
      obj@categories[[layer]] <- extra
    }
    changed <- c(changed, layer)
  }

  for (nm in unique(changed)) {
    v <- obj@layers[[nm]]
    if (!inherits(v, "rle")) {
      rv <- rle(v)
      if (utils::object.size(rv) < utils::object.size(v)) obj@layers[[nm]] <- rv
    }
  }
  if ("categories" %in% names(slots)) .check_categories(obj)

  if (!is.null(step) && step$record) {
    obj@provenance <- c(obj@provenance, list(.entry(step, activity)))
  }

  methods::validObject(obj)
  obj
}


#' Build the provenance entry of a step
#'
#' The entry follows the W3C PROV ontology: the activity
#' (\code{wasGeneratedBy}, with the call's arguments under
#' \code{withArguments}), the entity it read (\code{wasDerivedFrom}), the one
#' it produced (\code{generated}), the agent (\code{wasAssociatedWith}, the
#' package and version), the time and a digest of the activity. Inputs too
#' large to copy are listed under \code{used} with their digest. It is the same
#' shape the \pkg{bitfield} package writes.
#' @noRd

.entry <- function(step, activity = list()) {

  generated <- c(list(withArguments = step$args), activity)

  entry <- list(wasGeneratedBy    = generated,
                wasDerivedFrom    = if (is.null(step$from)) NA_character_ else step$from,
                generated         = if (is.null(step$to)) NA_character_ else step$to,
                wasAssociatedWith = step$agent,
                atTime            = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                hash              = digest::digest(generated))
  if (length(step$used)) entry$used <- step$used

  structure(list(entry), names = step$fn)
}


# arguments larger than this (bytes) are recorded as a digest, not copied
.max_arg_size <- 10000

# arguments that must NOT be evaluated when a step is recorded, because they
# are expressions over the layers of the mosaik the step will run on
.quoted_args <- list(mdf_filter = "expr",
                     mdf_loop   = "until")


# the name of the function a call invoked: the call's own name where it is one,
# otherwise found by looking the function up in its namespace (do.call with a
# function object, lapply)
.fn_name <- function(head, fun, ns) {
  if (is.call(head) && as.character(head[[1]]) %in% c("::", ":::")) {
    head <- head[[3]]
  }
  if (is.name(head) && identical(get0(as.character(head), envir = ns), fun)) {
    return(as.character(head))
  }
  for (nm in ls(ns, all.names = TRUE)) {
    if (identical(get(nm, envir = ns), fun)) return(nm)
  }
  deparse1(head)
}


# the software that ran a step: the package whose namespace holds the
# function, so a downstream package is recorded as itself
.agent <- function(ns) {
  if (!isNamespace(ns)) return(NA_character_)
  pkg <- getNamespaceName(ns)
  paste0(pkg, " ", utils::packageVersion(pkg))
}


# warn where a category table and the values of its layer disagree
.check_categories <- function(obj) {
  for (nm in names(obj@categories)) {
    entry <- obj@categories[[nm]]
    if (!"gid" %in% names(entry)) next
    v <- obj@layers[[nm]]
    if (inherits(v, "rle")) v <- inverse.rle(v)
    u <- sort(unique(v[!is.na(v)]))
    missing <- setdiff(u, entry$gid)
    if (length(missing)) {
      warning("Layer '", nm, "' contains cell values (",
              paste(missing, collapse = ", "), ") not listed in categories$gid.")
    }
    extra <- setdiff(entry$gid, u)
    if (length(extra)) {
      warning("Categories entry '", nm, "' has gid values (",
              paste(extra, collapse = ", "), ") not found in layer cell values.")
    }
  }
}


#' Read the arguments out of a provenance entry
#'
#' The one place that knows where arguments sit inside an entry, so that the
#' consumers (\code{show}, \code{\link{msk_vis}}, \code{\link{mdf}}) do not each
#' hard-code the path.
#' @noRd

.prov_args <- function(entry) entry$wasGeneratedBy$withArguments


#' Format a history for printing
#'
#' One line per step: its number, the function, the layer it read and the one
#' it wrote, and the arguments that differ from the function's defaults. The
#' layer arguments are left out of that list because the arrow already shows
#' them. Timestamps, the agent and the digests stay in the object.
#'
#' @param prov the \code{@provenance} of a mosaik.
#' @return character, one element per step.
#' @noRd

.format_history <- function(prov) {

  if (length(prov) == 0) return(character())

  rows <- lapply(prov, function(step) {
    if (is.character(step)) return(c(step, "", ""))
    fn <- names(step)[1]
    e <- step[[1]]

    from <- e$wasDerivedFrom
    to <- e$generated
    from <- if (is.null(from) || all(is.na(from))) "" else paste(from, collapse = ", ")
    to <- if (is.null(to) || all(is.na(to))) "" else to
    io <- if (nzchar(to) && to != from) paste0(from, if (nzchar(from)) " ", "-> ", to) else from

    args <- .prov_args(e)
    defaults <- .formals_of(fn, e$wasAssociatedWith)
    shown <- character()
    for (i in seq_along(args)) {
      nm <- names(args)[i]
      if (is.null(nm)) nm <- ""
      if (nm %in% c("layer", "layers", "add", "name")) next
      v <- args[[i]]
      if (nzchar(nm) && nm %in% names(defaults) &&
          .same_as_default(v, defaults[[nm]])) next
      shown <- c(shown, if (nzchar(nm)) paste0(nm, " = ", .format_arg(v)) else .format_arg(v))
    }
    if (!is.null(e$wasGeneratedBy$iterations)) {
      shown <- c(shown, paste0("(", e$wasGeneratedBy$iterations, " iterations)"))
    }
    c(fn, io, paste(shown, collapse = ", "))
  })

  rows <- do.call(rbind, rows)
  num <- formatC(seq_len(nrow(rows)), width = nchar(nrow(rows)))
  fn <- formatC(rows[, 1], width = -max(nchar(rows[, 1])))
  io <- formatC(rows[, 2], width = -max(nchar(rows[, 2])))
  trimws(paste(num, fn, io, rows[, 3], sep = "  "), which = "right")
}

# the formals of the function that wrote a step, from the package the entry
# names as its agent, so a mundus step is read against mundus
.formals_of <- function(fn, agent) {
  pkg <- if (is.character(agent) && !is.na(agent)) sub(" .*$", "", agent) else "mosaik"
  if (!requireNamespace(pkg, quietly = TRUE)) return(list())
  f <- get0(fn, envir = asNamespace(pkg), mode = "function")
  if (is.null(f)) list() else as.list(formals(f))
}

.same_as_default <- function(v, default) {
  if (is.null(v)) return(TRUE)
  d <- tryCatch(eval(default, baseenv()), error = function(e) NULL)
  identical(v, d)
}

# a short, readable form of one recorded argument
.format_arg <- function(v) {
  if (inherits(v, "msk_used")) return(paste0("<data ", substr(v$digest, 1, 8), ">"))
  if (methods::is(v, "struct")) return(paste0("<struct ", paste(dim(v@pattern), collapse = "x"), ">"))
  if (methods::is(v, "mosaik")) {
    n <- length(v@provenance)
    return(paste0("<recipe, ", n, if (n == 1) " step>" else " steps>"))
  }
  if (is.call(v) || is.name(v)) return(deparse1(v))
  if (is.function(v)) return("<function>")
  if (is.character(v) && length(v) == 1) return(paste0('"', v, '"'))
  if (is.atomic(v) && length(v) == 1) return(format(v))
  if (is.atomic(v) && length(v) <= 4) return(deparse1(v))
  if (is.atomic(v)) return(paste0("<", class(v)[1], " [", length(v), "]>"))
  paste0("<", class(v)[1], ">")
}
