# The default theme settings, and the set of valid theme keys.
#
# msk_vis() fills in from here whatever a .theme() call leaves unspecified, and
# .assert_theme_names() validates user settings against these names. Keeping
# the settings flat and dotted rather than nested is what lets one name check
# catch a typo at any depth, and what makes .theme()'s arguments completable.
#
# Adding a setting here makes it valid in .theme() automatically; document it
# in .theme()'s roxygen block as well, since that is where users read them.

.themeDefaults <- list(
  # titles ------------------------------------------------------------
  # figure.* heads the whole figure, panel.* heads each panel within it
  figure.title.plot     = TRUE,
  figure.title.fontsize = 14,
  figure.title.colour   = "black",
  panel.title.plot      = TRUE,
  panel.title.fontsize  = 12,
  panel.title.colour    = "black",

  # bounding box ------------------------------------------------------
  box.plot              = TRUE,
  box.linewidth         = 1,
  box.linetype          = "solid",
  box.linecol           = "black",

  # x axis ------------------------------------------------------------
  xAxis.plot            = TRUE,
  xAxis.bins            = 5,
  xAxis.margin          = 0.01,
  xAxis.label.plot      = TRUE,
  xAxis.label.title     = "x",
  xAxis.label.fontsize  = 10,
  xAxis.label.colour    = "grey45",
  xAxis.label.rotation  = 0,
  xAxis.ticks.plot      = TRUE,
  xAxis.ticks.fontsize  = 8,
  xAxis.ticks.colour    = "grey45",
  xAxis.ticks.rotation  = 0,
  xAxis.ticks.digits    = 1,

  # y axis ------------------------------------------------------------
  yAxis.plot            = TRUE,
  yAxis.bins            = 5,
  yAxis.margin          = 0.01,
  yAxis.label.plot      = TRUE,
  yAxis.label.title     = "y",
  yAxis.label.fontsize  = 10,
  yAxis.label.colour    = "grey45",
  yAxis.label.rotation  = 0,
  yAxis.ticks.plot      = TRUE,
  yAxis.ticks.fontsize  = 8,
  yAxis.ticks.colour    = "grey45",
  yAxis.ticks.rotation  = 0,
  yAxis.ticks.digits    = 1,

  # grid --------------------------------------------------------------
  grid.plot             = TRUE,
  grid.colour           = "grey70",
  grid.linetype         = "dashed",
  grid.linewidth        = 0.5,

  # legend ------------------------------------------------------------
  legend.plot           = TRUE,
  legend.bins           = 5,
  legend.ascending      = TRUE,
  legend.position       = "right",
  legend.yRatio         = 0.7,
  legend.digits         = 1,
  legend.label.plot     = TRUE,
  legend.label.fontsize = 8,
  legend.label.colour   = "grey45",
  legend.box.plot       = TRUE,
  legend.box.linetype   = "solid",
  legend.box.linewidth  = 1,
  legend.box.colour     = "black",

  # colour scale ------------------------------------------------------
  # colours takes ramp stops, a palette name (see .resolve_colours), or a
  # named vector mapping category labels to colours.
  colours               = c("#00204D", "#FFEA46"),
  colours.missing       = "#FFFFFF00",
  scale.range           = NULL,
  scale.bins            = 256
)


#' Validate theme keys against the known paths
#'
#' The theme is a flat set of dotted paths, so a typo at any depth is caught by
#' one name check -- there is no nesting to walk. Unknown keys are an error, not
#' a warning: a silently ignored setting looks like the plot simply does not
#' respond to it, which is far harder to diagnose than a failed call.
#'
#' @param x [named list][list]\cr the user-supplied theme.
#' @param known [character][character]\cr the valid paths.
#' @return invisibly \code{TRUE}; errors when a key is unknown.
#' @importFrom utils adist
#' @noRd

.assert_theme_names <- function(x, known = names(.themeDefaults)){

  if(is.null(x)) return(invisible(TRUE))

  nms <- names(x)
  if(length(x) > 0 && (is.null(nms) || any(nms == ""))){
    stop("every theme setting must be named, e.g. .theme(legend.plot = FALSE).",
         call. = FALSE)
  }

  bad <- setdiff(nms, known)
  if(length(bad) == 0) return(invisible(TRUE))

  # suggest the closest known key: most unknown names are typos or a guess at
  # a name that follows a different convention, and both are fixed by seeing
  # the intended spelling.
  hints <- vapply(bad, function(b){
    d <- adist(b, known, ignore.case = TRUE)[1, ]
    if(min(d) <= max(3, nchar(b) %/% 3)){
      paste0("'", b, "' -- did you mean '", known[which.min(d)], "'?")
    } else {
      paste0("'", b, "'")
    }
  }, character(1))

  stop("unknown theme setting(s):\n  ", paste(hints, collapse = "\n  "),
       "\nsee ?.theme for the full list of settings.", call. = FALSE)
}


#' Set the graphical parameters of a plot
#'
#' Collects the settings that govern how a \code{\link{msk_vis}} plot is drawn:
#' its titles, axes, grid, legend and default colour scale. Pass the result to
#' the \code{theme} argument.
#'
#' Each setting is a single dotted name, so they can be typed straight into the
#' call and completed by name:
#'
#' \preformatted{
#'   msk_vis(obj, theme = .theme(legend.plot = FALSE))
#'
#'   mine <- .theme(colours = "viridis", panel.title.fontsize = 14)
#'   msk_vis(obj, theme = mine)
#'   msk_vis(obj, theme = .theme(mine, legend.plot = FALSE))
#' }
#'
#' Anything not given keeps its default. Settings are:
#'
#' \strong{titles} -- \code{figure.title.plot}, \code{figure.title.fontsize},
#'   \code{figure.title.colour} for the heading over the whole figure (set its
#'   text with \code{msk_vis}'s \code{title} argument), and
#'   \code{panel.title.plot}, \code{panel.title.fontsize},
#'   \code{panel.title.colour} for the heading over each panel (set its text
#'   with \code{\link{.layer}}'s \code{panel} argument).
#'
#' \strong{box} (the frame around the raster) -- \code{box.plot},
#'   \code{box.linewidth}, \code{box.linetype}, \code{box.linecol}.
#'
#' \strong{axes} -- for both \code{xAxis} and \code{yAxis}: \code{.plot},
#'   \code{.bins} (number of tick intervals), \code{.margin}, and the
#'   sub-groups \code{.label.plot}, \code{.label.title}, \code{.label.fontsize},
#'   \code{.label.colour}, \code{.label.rotation}, \code{.ticks.plot},
#'   \code{.ticks.fontsize}, \code{.ticks.colour}, \code{.ticks.rotation},
#'   \code{.ticks.digits}.
#'
#' \strong{grid} -- \code{grid.plot}, \code{grid.colour}, \code{grid.linetype},
#'   \code{grid.linewidth}.
#'
#' \strong{legend} -- \code{legend.plot}, \code{legend.bins} (number of labelled
#'   ticks), \code{legend.ascending}, \code{legend.position},
#'   \code{legend.yRatio}, \code{legend.digits}, \code{legend.label.plot},
#'   \code{legend.label.fontsize}, \code{legend.label.colour},
#'   \code{legend.box.plot}, \code{legend.box.linetype},
#'   \code{legend.box.linewidth}, \code{legend.box.colour}.
#'
#' \strong{colour scale} -- \code{colours} (ramp stops, a palette name, or a
#'   named vector mapping categories to colours), \code{colours.missing} (the
#'   colour for NA cells), \code{scale.range} (forced limits) and
#'   \code{scale.bins} (number of colour steps). A \code{\link{.layer}} can
#'   override any of these for itself.
#'
#' @param ... named settings from the list above. An unknown name is an error,
#'   with the nearest valid name suggested.
#' @param from \code{mskTheme}\cr an existing theme to extend, given as the
#'   first argument. Settings in \code{...} replace its values.
#' @return a list of class \code{mskTheme}.
#' @examples
#' .theme(legend.plot = FALSE)
#' .theme(colours = "viridis", panel.title.fontsize = 14)
#' @export

.theme <- function(from = NULL, ...){

  settings <- list(...)
  .assert_theme_names(settings)

  if(!is.null(from) && !inherits(from, "mskTheme")){
    stop("'from' must be a theme built with .theme().", call. = FALSE)
  }

  out <- if(is.null(from)) list() else unclass(from)
  out[names(settings)] <- settings

  structure(out, class = c("mskTheme", "list"))
}


#' Complete a user theme with the defaults
#'
#' @param theme \code{mskTheme}\cr user settings, or \code{NULL}.
#' @return the full theme list.
#' @noRd

.complete_theme <- function(theme = NULL){

  out <- .themeDefaults
  if(is.null(theme)) return(out)

  if(!inherits(theme, "mskTheme")){
    stop("'theme' must be built with .theme(), e.g. ",
         "theme = .theme(legend.plot = FALSE).", call. = FALSE)
  }

  theme <- unclass(theme)
  .assert_theme_names(theme)
  out[names(theme)] <- theme
  out
}


#' Resolve a colour specification into a ramp or a category map
#'
#' Three forms, told apart without a \code{type} argument:
#' \itemize{
#'   \item a \strong{named} character vector maps category labels to colours,
#'     e.g. \code{c(forest = "darkgreen", crop = "khaki")};
#'   \item a \strong{single string} that names a palette is passed to
#'     \code{grDevices::hcl.colors}, e.g. \code{"viridis"} or \code{"terrain"};
#'   \item any \strong{other} character vector is taken as ramp stops and
#'     interpolated, e.g. \code{c("black", "white")}.
#' }
#'
#' @param colours the \code{colours} theme setting.
#' @param n [numeric(1)][numeric]\cr how many colours the ramp needs.
#' @return a character vector of \code{n} colours, or -- for the named form --
#'   the named vector unchanged, for the caller to match against categories.
#' @importFrom grDevices hcl.colors hcl.pals colorRampPalette
#' @noRd

.resolve_colours <- function(colours, n){

  # named vector: a category map, returned as-is
  if(!is.null(names(colours))) return(colours)

  if(length(colours) == 1 && is.character(colours)){
    # hcl.pals() capitalises ("Viridis", "Terrain"); match case-insensitively so
    # the spelling a user naturally reaches for works
    hit <- match(tolower(colours), tolower(hcl.pals()))
    if(!is.na(hit)){
      return(hcl.colors(n = n, palette = hcl.pals()[hit]))
    }
    # a single unnamed colour is a legitimate one-stop "ramp"
    if(!colours %in% grDevices::colors() &&
       !grepl("^#[0-9A-Fa-f]{6,8}$", colours)){
      stop("'", colours, "' is neither a colour nor a known palette. ",
           "See grDevices::hcl.pals() for the available palettes.",
           call. = FALSE)
    }
  }

  colorRampPalette(colors = colours)(n)
}
