#' Specify one layer of a plot
#'
#' Describes how a single layer is drawn: which colours it uses, whether it is
#' shaded by relief, and which panel it belongs to. Pass as many of these to
#' \code{\link{msk_vis}} as the plot needs -- one per layer, in drawing order.
#'
#' Layers sharing a \code{panel} are composited into that one panel, bottom to
#' top in the order given, with NA cells transparent. This is how a variable is
#' drawn over terrain: give the terrain and the variable the same \code{panel},
#' and each keeps its own colour scale and its own legend.
#'
#' @param layer [`character(1)`][character]\cr name of the layer to draw.
#' @param panel [`character(1)`][character]\cr which panel this layer belongs to.
#'   Layers sharing a value are drawn into the same panel, and the value is that
#'   panel's title. \code{NULL} (default) gives the layer a panel of its own,
#'   titled after the layer.
#' @param colours [`character`][character]\cr the layer's colour scale: ramp stops
#'   (\code{c("black", "white")}), a palette name
#'   (\code{"viridis"}; see \code{grDevices::hcl.pals}), or a \strong{named}
#'   vector mapping category labels to colours
#'   (\code{c(forest = "darkgreen", crop = "khaki")}). \code{NULL} takes the
#'   colours of the layer's categories where they carry a \code{colour} field,
#'   else the theme's \code{colours}.
#' @param limits [`numeric(2)`][numeric]\cr \code{c(min, max)} forcing this
#'   layer's colour scale, instead of taking the range from its own values.
#'   Cells outside are clamped to the end colours rather than dropped.
#' @param bins [`numeric(1)`][numeric]\cr number of colour steps for a continuous
#'   layer. \code{NULL} takes the theme's \code{scale.bins}.
#' @param hillshade [`named list`][list]\cr shade this layer by terrain relief.
#'   The shading multiplies the layer's own colours, so any layer can be draped
#'   over the relief. \code{NULL} (default) disables it. Elements:
#'   \itemize{
#'     \item \code{layer} [`character(1)`][character] name of the elevation layer
#'       to shade by (required).
#'     \item \code{azimuth} [`numeric`][numeric] light direction in degrees
#'       clockwise from north, default 315. Give several to average them into a
#'       multi-directional hillshade, e.g. \code{c(225, 270, 315, 360)}, which
#'       keeps slopes facing away from a single light from going flat.
#'     \item \code{altitude} [`numeric(1)`][numeric] light height in degrees above
#'       the horizon, default 45.
#'     \item \code{exaggeration} [`numeric(1)`][numeric] vertical exaggeration of
#'       the elevation before the gradient is taken, default 1.
#'     \item \code{intensity} [`numeric(1)`][numeric] strength of the effect in
#'       \code{[0, 1]}, default 1.
#'   }
#' @param legend [`logical(1)`][logical]\cr give this layer its own legend
#'   (\code{TRUE}, default). Set \code{FALSE} for a layer that needs no scale of
#'   its own, such as a hillshade base.
#' @param title [`character(1)`][character]\cr heading for this layer's legend.
#'   \code{NULL} uses the layer name.
#' @return a list of class \code{mskLayer}.
#' @seealso \code{\link{msk_vis}}, which draws these; \code{\link{.theme}} for
#'   settings that apply to the whole plot.
#' @examples
#' # the plainest form: just a layer name
#' msk_vis(landscape, .layer("cover"))
#'
#' # a palette, and a scale fixed so the figure is comparable across calls
#' msk_vis(landscape, .layer("cover", colours = "viridis", limits = c(0, 60)))
#'
#' # an explicit colour per category, which fixes the meaning of each colour
#' # across every figure in a series
#' m <- mdf_filter(landscape, cover == 47, add = "forest")
#' msk_vis(m, .layer("forest", colours = c("0" = "khaki", "1" = "darkgreen")))
#'
#' # a variable drawn over shaded terrain: same panel, one legend each, and
#' # the flagged cells stay saturated because only the terrain carries the
#' # shading
#' t <- mosaik(extent = c(0, 10000, 0, 10000), res = 100) |>
#'   drw_texture(name = "dem", hurst = 0.9, seed = 1) |>
#'   mdf_scale(layer = "dem", range = c(0, 200))
#' t <- mdf_filter(t, dem >= 160, add = "high")
#' t <- mdf_filter(t, high == 1, value = TRUE, layer = "high")
#' msk_vis(t,
#'         .layer("dem", panel = "relief", colours = "terrain",
#'                hillshade = list(layer = "dem", exaggeration = 5)),
#'         .layer("high", panel = "relief", colours = "red", legend = FALSE))
#' @importFrom checkmate assertCharacter assertNumeric assertNumber assertList
#'   assertSubset assertLogical
#' @export

.layer <- function(layer, panel = NULL, colours = NULL, limits = NULL,
                   bins = NULL, hillshade = NULL, legend = TRUE,
                   title = NULL){

  assertCharacter(x = layer, len = 1, any.missing = FALSE)
  assertCharacter(x = colours, min.len = 1, any.missing = FALSE, null.ok = TRUE)
  assertNumeric(x = limits, len = 2, any.missing = FALSE, sorted = TRUE,
                finite = TRUE, null.ok = TRUE)
  assertNumber(x = bins, lower = 1, null.ok = TRUE)
  assertLogical(x = legend, len = 1, any.missing = FALSE)
  assertCharacter(x = title, len = 1, any.missing = FALSE, null.ok = TRUE)
  assertCharacter(x = panel, len = 1, any.missing = FALSE, null.ok = TRUE)
  if(!is.null(hillshade)){
    assertList(x = hillshade, names = "unique")
    assertSubset(x = names(hillshade),
                 choices = c("layer", "azimuth", "altitude", "exaggeration",
                             "intensity"))
    if(is.null(hillshade$layer)){
      stop("hillshade needs a 'layer' element naming the elevation layer to ",
           "shade by, e.g. hillshade = list(layer = \"terrain\").", call. = FALSE)
    }
    assertCharacter(x = hillshade$layer, len = 1, any.missing = FALSE)
    if(!is.null(hillshade$intensity)){
      assertNumber(x = hillshade$intensity, lower = 0, upper = 1)
    }
  }

  structure(list(layer = layer, panel = panel, colours = colours,
                 limits = limits, bins = bins, hillshade = hillshade,
                 legend = legend, title = title),
            class = c("mskLayer", "list"))
}


#' Visualise a mosaik
#'
#' Draws the layers of one mosaik. What is drawn is described by
#' \code{\link{.layer}} specifications, and how it is drawn by a
#' \code{\link{.theme}}.
#'
#' @details
#' \strong{Panels and compositing.} Each \code{.layer()} gets a panel of its
#' own unless it is given a \code{panel}; layers sharing a \code{panel} value
#' are drawn into that one panel, in the order they appear, with NA cells
#' transparent. This is the difference between showing two variables side by
#' side and drawing one on top of the other:
#'
#' \preformatted{
#'   # two panels
#'   msk_vis(m, .layer("dem"), .layer("rivers"))
#'
#'   # one panel, rivers over terrain
#'   msk_vis(m, .layer("dem",    panel = "catchment"),
#'              .layer("rivers", panel = "catchment"))
#' }
#'
#' Every layer keeps its own colour scale and gets its own legend, so a variable
#' drawn over terrain is read against its own numbers rather than the terrain's.
#' The \code{panel} value is also the panel's title; \code{title} heads the
#' figure as a whole.
#'
#' Only NA cells let the layer beneath show through, so a layer with a value
#' everywhere hides whatever it is drawn on. An overlay is therefore a layer
#' that is mostly NA -- rivers on a bare grid, or a subset masked out of a full
#' layer with \code{\link{mdf_filter}}. Set \code{legend = FALSE} on a layer that
#' flags cells rather than carrying a scale of its own.
#'
#' \strong{Comparing panels.} By default each panel scales its colours to its
#' own values, which means an identical colour in two panels can stand for
#' different numbers. \code{shared_scale = TRUE} pools the range per layer name
#' so the panels become comparable; \code{limits} on a \code{.layer()} fixes a
#' range outright, which is what makes a figure comparable across separate
#' calls.
#'
#' \strong{Relief.} A layer can be shaded by terrain via the \code{hillshade}
#' element of \code{.layer()}. The shading multiplies that layer's own colours,
#' so any variable can be draped over the relief, and a layer composited on top
#' of a shaded one is not darkened by it.
#'
#' @param obj [`mosaik`][mosaik]\cr the object to plot.
#' @param ... \code{mskLayer}\cr layer specifications built with
#'   \code{\link{.layer}}, one per layer to draw. With none given, every layer
#'   of \code{obj} is drawn in a panel of its own.
#' @param title [`character(1)`][character]\cr heading for the figure as a whole,
#'   drawn above the panels. \code{NULL} (default) draws none and gives the
#'   panels the full page. Panel headings come from \code{\link{.layer}}'s
#'   \code{panel} argument instead.
#' @param shared_scale [`logical(1)`][logical]\cr give every panel showing the
#'   same layer one common colour scale, spanning that layer's range across all
#'   panels (\code{TRUE}). Layers with an explicit \code{limits} keep it.
#'   Default \code{FALSE}.
#' @param window [`numeric(4)`][numeric]\cr extent as \code{c(xmin, xmax, ymin,
#'   ymax)} to which the plot is limited. If \code{NULL}, the object's extent is
#'   used.
#' @param theme \code{mskTheme}\cr graphical settings built with
#'   \code{\link{.theme}}, e.g. \code{.theme(legend.plot = FALSE,
#'   panel.title.fontsize = 14)}. Anything left out keeps its default.
#' @param trace [`logical(1)`][logical]\cr Print the provenance information of the
#'   mosaik object (\code{TRUE}), or simply plot (\code{FALSE}, default).
#' @return Returns invisibly an object of class \code{recordedplot}.
#' @seealso \code{\link{.layer}} for what to draw, \code{\link{.theme}} for how
#'   to draw it.
#' @examples
#' # every layer, one panel each
#' msk_vis(landscape)
#'
#' # one named layer, with a palette of its own
#' msk_vis(landscape, .layer("cover", colours = "viridis"))
#'
#' # two layers on a common scale, so the panels can be compared
#' msk_vis(landscape, .layer("cover"), .layer("canopy"),
#'         shared_scale = TRUE)
#'
#' # two layers in one panel: an overlay only shows through where the upper
#' # layer is NA, so it is built by masking down to the cells of interest --
#' # here the cells without vegetation, flagged over the cover map
#' m <- mdf_filter(landscape, canopy == 0, add = "flat")
#' m <- mdf_filter(m, flat == 1, value = TRUE, layer = "flat")
#' msk_vis(m,
#'         .layer("cover",  panel = "landscape", colours = "Terrain 2"),
#'         .layer("flat", panel = "landscape", colours = "red",
#'                legend = FALSE))
#'
#' # a heading over the whole figure, panels headed individually beneath it
#' msk_vis(landscape, .layer("cover", panel = "what grows"),
#'         .layer("canopy", panel = "how tall it grows"),
#'         title = "land cover and vegetation")
#'
#' # drop the chrome
#' msk_vis(landscape, theme = .theme(panel.title.plot = FALSE,
#'                                   grid.plot = FALSE))
#' @importFrom checkmate assertClass assertLogical assertNumeric assertCharacter
#'   assertList assertSubset assertNumber
#' @importFrom grid grid.newpage pushViewport viewport grid.layout grid.rect
#'   grid.raster grid.clip unit grid.draw grid.grill upViewport grid.text gpar
#'   unit.c convertX
#' @importFrom grDevices recordPlot dev.list
#' @export

msk_vis <- function(obj, ..., title = NULL, shared_scale = FALSE,
                     window = NULL, theme = NULL, trace = FALSE){

  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = title, len = 1, any.missing = FALSE, null.ok = TRUE)
  assertLogical(x = shared_scale, len = 1, any.missing = FALSE)
  assertLogical(x = trace, len = 1, any.missing = FALSE)
  theme <- .complete_theme(theme)

  specs <- list(...)
  if(length(specs) > 0){
    isSpec <- vapply(specs, inherits, logical(1), what = "mskLayer")
    if(!all(isSpec)){
      stop("every argument after 'obj' must be a layer specification built ",
           "with .layer(); element(s) ", paste(which(!isSpec), collapse = ", "),
           " are not.", call. = FALSE)
    }
  } else {
    # no specs: every layer in a panel of its own, the plain overview
    specs <- lapply(names(obj@layers), .layer)
  }

  if(length(specs) == 0){
    stop("'obj' has no layers to plot.", call. = FALSE)
  }

  missingLayers <- setdiff(vapply(specs, function(s) s$layer, character(1)),
                           names(obj@layers))
  if(length(missingLayers) > 0){
    stop("layer(s) not in this mosaik: ",
         paste0("'", missingLayers, "'", collapse = ", "), ".\navailable: ",
         paste0("'", names(obj@layers), "'", collapse = ", "), call. = FALSE)
  }

  # group specs into panels: an explicit panel value groups, NULL stands alone.
  # Using the layer name as the implicit key would merge two specs of the same
  # layer, which is a legitimate thing to want in separate panels.
  panelKey <- vapply(seq_along(specs), function(i){
    p <- specs[[i]]$panel
    if(is.null(p)) paste0("\r__solo", i) else p
  }, character(1))
  panelOrder <- unique(panelKey)

  # a panel's title is its key when that key was named, else the layer names
  panelTitles <- vapply(panelOrder, function(k){
    members <- specs[panelKey == k]
    if(!grepl("^\r__solo", k)){
      k
    } else {
      paste0(vapply(members, function(s) s$layer, character(1)), collapse = " / ")
    }
  }, character(1))

  panels <- length(panelOrder)

  if(panels > 15){
    message(paste0("this will produce ", panels, " panels."))
  }

  # panel grid

  # up to three panels side by side, more in a square grid
  if(panels <= 3){
    ncol <- panels
  } else {
    ncol <- ceiling(sqrt(panels))
  }
  nrow <- ceiling(panels / ncol)
  panelPosY <- rep(rev(seq(from = 1, to = nrow)), each = ncol)
  panelPosX <- rep(seq(from = 1, to = ncol), times = nrow)

  # shared scale: pool each layer's values across every panel that shows it, so
  # that one colour means one value throughout. Pooling is per layer name --
  # putting elevation and temperature on a common scale would be meaningless.
  sharedLimits <- list()
  if(shared_scale){
    for(lyr in unique(vapply(specs, function(s) s$layer, character(1)))){
      vals <- msk_pull(obj, lyr)
      vals <- vals[!is.na(vals)]
      if(!is.numeric(vals) || length(vals) == 0) next
      rng <- range(vals)
      # a layer that is constant or non-numeric keeps its per-panel scale
      if(is.finite(rng[1]) && rng[2] > rng[1]) sharedLimits[[lyr]] <- rng
    }
  }

  grid.newpage()
  pushViewport(viewport(name = "mosaik"))

  # a figure title takes a band off the top of the page and the panels share
  # what is left; with no title nothing is reserved, so the panel grid fills
  # the page exactly as it does without one
  if(!is.null(title) && theme$figure.title.plot){
    bandH <- unit(theme$figure.title.fontsize + 10, "points")
    grid.text(label = title,
              y = unit(1, "npc") - unit(5, "points"),
              just = "top",
              gp = gpar(fontsize = theme$figure.title.fontsize,
                        col = theme$figure.title.colour),
              name = "figureTitle")
    pushViewport(viewport(y = 0, just = "bottom",
                          height = unit(1, "npc") - bandH,
                          name = "panels"))
  }

  # determine window once: all panels come from one mosaik, so they share it
  # the cells are drawn over the whole window, so a window smaller than the
  # map needs the map cropped to it first, snapped outward to whole cells
  if(!is.null(window)){
    assertNumeric(x = window, len = 4)
    # only the cells are drawn, so the measured values go before cropping
    obj@categories <- lapply(obj@categories, function(e) e[setdiff(names(e), .measured_names(e))])
    obj@global <- list()
    obj <- mdf_crop(obj, extent = window)
  }
  ext <- obj@extent
  theWindow <- list(xmin = ext[1], xmax = ext[2],
                    ymin = ext[3], ymax = ext[4])

  for(i in seq_along(panelOrder)){

    theName <- panelTitles[i]
    theSpecs <- specs[panelKey == panelOrder[i]]

    temp <- .makePlot(obj = obj, specs = theSpecs,
                      sharedLimits = sharedLimits,
                      window = theWindow, theme = theme)
    theTheme <- temp$theme
    theGrob <- temp$grob
    theLegend <- temp$legend
    theLayout <- temp$layout

    # create panel viewport
    pushViewport(viewport(x = (panelPosX[i] / ncol) - (1 / ncol / 2),
                          y = (panelPosY[i] / nrow) - (1 / nrow / 2),
                          width = 1 / ncol,
                          height = 1 / nrow,
                          name = theName))
    grid.rect(width = convertX(unit(1, "npc"), "native"),
              gp = gpar(col = NA, fill = NA), name = "panelGrob")

    grid.rect(x = unit(theLayout$window$xmin, "points"),
              y = unit(theLayout$window$ymin, "points"),
              height = unit(theLayout$window$ymax - theLayout$window$ymin, "points"),
              width = unit(theLayout$window$xmax - theLayout$window$xmin, "points"),
              gp = gpar(col = NA, fill = NA), name = "extentGrob")

    # plot layout
    myLayout <- grid.layout(nrow = 4, ncol = 3,
                            widths = unit.c(unit(theLayout$dim$x1, "points"),
                                            theLayout$dim$x2,
                                            unit(theLayout$dim$x3, "points")),
                            heights = unit.c(unit(theLayout$dim$y1, "points"),
                                             theLayout$dim$y2,
                                             unit(theLayout$dim$y3, "points"),
                                             unit(theLayout$dim$y4, "points")))

    layoutVP <- viewport(name = "theLayout", layout = myLayout)
    titleVP <- viewport(name = "title",
                        layout.pos.col = 2, layout.pos.row = 1)
    yAxisVP <- viewport(name = "y_axis",
                        layout.pos.col = 1, layout.pos.row = 2,
                        xscale = c(theLayout$scale$xmin, theLayout$scale$xmax),
                        yscale = c(theLayout$scale$ymin, theLayout$scale$ymax))
    xAxisVP <- viewport(name = "x_axis",
                        layout.pos.col = 2, layout.pos.row = 3,
                        xscale = c(theLayout$scale$xmin, theLayout$scale$xmax),
                        yscale = c(theLayout$scale$ymin, theLayout$scale$ymax))
    legendVP <- viewport(name = "legend",
                         layout.pos.col = theLayout$legend$posX,
                         layout.pos.row = theLayout$legend$posY)
    plotVP <- viewport(name = "grid",
                       layout.pos.col = 2, layout.pos.row = 2,
                       xscale = c(theLayout$scale$xmin, theLayout$scale$xmax),
                       yscale = c(theLayout$scale$ymin, theLayout$scale$ymax))
    boxVP <- viewport(width = unit(1, "npc") - unit(2 * theLayout$margin$x, "native") + unit(theTheme$box.linewidth, "points"),
                      height = unit(1, "npc") - unit(2 * theLayout$margin$y, "native") + unit(theTheme$box.linewidth, "points"),
                      xscale = c(theLayout$scale$xmin, theLayout$scale$xmax),
                      yscale = c(theLayout$scale$ymin, theLayout$scale$ymax),
                      name = "box")
    legBoxVP <- viewport(height = unit(1, "npc") * theTheme$legend.yRatio,
                         width = unit(1, "npc"),
                         name = "legend_box")

    pushViewport(layoutVP)

    # panel title
    if(theTheme$panel.title.plot){
      pushViewport(titleVP)
      grid.text(y = unit(1, "npc") - unit(3, "points"),
                just = "top",
                label = theName,
                gp = gpar(fontsize = theTheme$panel.title.fontsize,
                          col = theTheme$panel.title.colour))
      upViewport()
    }

    # y axis
    if(theTheme$yAxis.plot){
      pushViewport(yAxisVP)
      if(theTheme$yAxis.label.plot){
        grid.text(x = unit(1, "npc") - unit(5, "points") - unit(theLayout$labels$yAxisTicksW, "points"),
                  just = "right",
                  label = theTheme$yAxis.label.title,
                  rot = theTheme$yAxis.label.rotation,
                  name = "y_title",
                  gp = gpar(fontsize = theTheme$yAxis.label.fontsize,
                            col = theTheme$yAxis.label.colour))
      }
      if(theTheme$yAxis.ticks.plot){
        grid.text(x = unit(1, "npc") - unit(2, "points"),
                  label = as.character(round(theLayout$grid$yMaj, theTheme$yAxis.ticks.digits)),
                  just = "right",
                  y = unit(theLayout$grid$yMaj, "native"),
                  rot = theTheme$xAxis.ticks.rotation,
                  name = "y_tick_labels",
                  gp = gpar(fontsize = theTheme$yAxis.ticks.fontsize,
                            col = theTheme$yAxis.ticks.colour))
      }
      upViewport()
    }

    # x axis
    if(theTheme$xAxis.plot){
      pushViewport(xAxisVP)
      if(theTheme$xAxis.label.plot){
        grid.text(y = unit(1, "npc") - unit(3, "points") - unit(theLayout$labels$xAxisTicksH, "points"),
                  just = "top",
                  label = theTheme$xAxis.label.title,
                  rot = theTheme$xAxis.label.rotation,
                  name = "x_title",
                  gp = gpar(fontsize = theTheme$xAxis.label.fontsize,
                            col = theTheme$xAxis.label.colour))
      }
      if(theTheme$xAxis.ticks.plot){
        grid.text(label = as.character(round(theLayout$grid$xMaj, theTheme$xAxis.ticks.digits)),
                  x = unit(theLayout$grid$xMaj, "native"),
                  y = unit(1, "npc") - unit(theLayout$labels$xAxisTicksH, "points"),
                  just = "bottom",
                  rot = theTheme$xAxis.ticks.rotation,
                  name = "x_tick_labels",
                  gp = gpar(fontsize = theTheme$xAxis.ticks.fontsize,
                            col = theTheme$xAxis.ticks.colour))
      }
      upViewport()
    }

    # legend
    if(theTheme$legend.plot){
      pushViewport(legendVP)
      pushViewport(legBoxVP)
      for(j in seq_along(theLegend)){
        grid.draw(theLegend[[j]])
      }
      upViewport()
      upViewport()
    }

    # plot area
    pushViewport(plotVP)

    if(theTheme$grid.plot){
      grid.grill(h = unit(theLayout$grid$yMaj, "native"),
                 v = unit(theLayout$grid$xMaj, "native"),
                 gp = gpar(col = theTheme$grid.colour,
                           lwd = theTheme$grid.linewidth,
                           lty = theTheme$grid.linetype))
    }

    pushViewport(boxVP)
    if(theTheme$box.plot){
      grid.clip(width = unit(1, "npc"),
                height = unit(1, "npc"))
      grid.draw(theGrob)
      upViewport() # exit boxVP
    }

    upViewport() # exit plotVP
    upViewport() # exit layoutVP
    upViewport() # exit panel VP
  }

  # provenance trace: one mosaik, so this is printed once, not once per panel
  if(trace){
    steps <- .format_history(obj@provenance)
    if(length(steps) > 0){
      message("this object has the following history:\n  ",
              paste0(steps, collapse = "\n  "))
    }
  }

  upViewport() # exit 'mosaik' VP

  invisible(recordPlot(attach = "mosaik"))
}

#' Get the number of decimal places
#'
#' @param x [`numeric(1)`][numeric]\cr the number for which to derive decimal
#'   places.
#' @noRd

.getDecimals <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


#' Make the layout of a plot
#'
#' @param obj a mosaik object.
#' @param legend the legend object built with \code{.makeLegend}.
#' @param window named list with \code{xmin, xmax, ymin, ymax}.
#' @param theme the completed theme list.
#' @importFrom grid convertX unit
#' @noRd

.makeLayout <- function(obj, legend, window, theme){

  if(is.null(theme)) theme <- .themeDefaults

  maxWinX <- window$xmax
  minWinX <- window$xmin
  maxWinY <- window$ymax
  minWinY <- window$ymin

  xBins <- theme$xAxis.bins
  yBins <- theme$yAxis.bins

  ratio <- list(x = (maxWinX - minWinX) / (maxWinY - minWinY),
                y = (maxWinY - minWinY) / (maxWinX - minWinX))
  xBinSize <- (maxWinX - minWinX) / xBins
  yBinSize <- (maxWinY - minWinY) / yBins
  axisSteps <- list(x1 = seq(from = minWinX, to = maxWinX,
                             by = (maxWinX - minWinX) / xBins),
                    x2 = seq(from = minWinX + (xBinSize / 2), to = maxWinX,
                             by = (maxWinX - minWinX) / xBins),
                    y1 = seq(from = minWinY, to = maxWinY,
                             by = (maxWinY - minWinY) / yBins),
                    y2 = seq(from = minWinY + (yBinSize / 2), to = maxWinY,
                             by = (maxWinY - minWinY) / yBins))
  margin <- list(x = (maxWinX - minWinX) * theme$yAxis.margin,
                 y = (maxWinY - minWinY) * theme$xAxis.margin)

  titleH <- if(theme$panel.title.plot) theme$panel.title.fontsize + 6 else 0

  if(theme$legend.plot){
    legendW <- 0
    for(i in seq_along(legend)){
      labels <- legend[[i]]$children$legend_labels$label
      maxLbl <- labels[which.max(nchar(labels))]
      tempW <- as.numeric(ceiling(convertX(unit(1, "strwidth", maxLbl) + unit(25, "points"), "points")))
      legendW <- legendW + tempW
    }
    legendW <- unit(legendW, "points")
    legendH <- unit(0, "points")
  } else {
    legendW <- unit(0, "points")
    legendH <- unit(0, "points")
  }

  if(theme$legend.position == "right"){
    legendPosX <- 3
    legendPosY <- 2
  } else {
    legendPosX <- 2
    legendPosY <- 4
  }

  if(theme$yAxis.plot){
    yAxisTitleW <- theme$yAxis.label.fontsize + 5
    digits <- round(axisSteps$y1, theme$yAxis.ticks.digits)
    yAxisTicksW <- ceiling(convertX(unit(1, "strwidth", as.character(digits[which.max(nchar(digits))])), "points"))
    yAxisTicksW <- as.numeric(yAxisTicksW)
  } else {
    yAxisTitleW <- 0
    yAxisTicksW <- 0
  }
  if(theme$xAxis.plot){
    xAxisTitleH <- theme$xAxis.label.fontsize + 2
    xAxisTicksH <- theme$xAxis.ticks.fontsize
  } else {
    xAxisTitleH <- 0
    xAxisTicksH <- 0
  }

  gridH <- unit(1, "grobheight", "panelGrob") - unit(xAxisTitleH, "points") - unit(xAxisTicksH, "points") - unit(titleH, "points")
  gridW <- unit(1, "grobwidth", "panelGrob") - unit(yAxisTitleW, "points") - unit(yAxisTicksW, "points") - unit(legendW, "points")
  gridHr <- gridW * ratio$y
  gridWr <- gridH * ratio$x
  gridH <- min(gridH, gridHr)
  gridW <- min(gridW, gridWr)

  list(dim = list(x1 = yAxisTitleW + yAxisTicksW, x2 = gridW, x3 = legendW,
                  y1 = titleH, y2 = gridH, y3 = xAxisTitleH + xAxisTicksH, y4 = legendH),
       margin = list(x = margin$x, y = margin$y),
       window = list(xmin = minWinX, xmax = maxWinX,
                     ymin = minWinY, ymax = maxWinY),
       labels = list(titleH = titleH,
                     yAxisTicksW = yAxisTicksW,
                     xAxisTicksH = xAxisTicksH),
       scale = list(xmin = minWinX - margin$x, xmax = maxWinX + margin$x,
                    ymin = minWinY - margin$y, ymax = maxWinY + margin$y),
       grid = list(xMaj = axisSteps$x1, xMin = axisSteps$x2,
                   yMaj = axisSteps$y1, yMin = axisSteps$y2),
       legend = list(posX = legendPosX, posY = legendPosY))
}


#' Make one legend
#'
#' Builds the legend for a single layer. A panel holding several layers calls
#' this once per layer and offsets each by \code{prevX}, so every scale in the
#' panel is documented by its own strip.
#'
#' @param scaleValues the scale values.
#' @param colours the resolved colour vector for this layer.
#' @param title [`character(1)`][character]\cr heading above the strip, or NULL.
#' @param prevX [`unit`][unit]\cr horizontal offset, so successive legends in one
#'   panel sit side by side rather than on top of each other.
#' @param theme the completed theme list.
#' @importFrom grid textGrob rasterGrob rectGrob gpar gTree gList unit convertX
#' @importFrom grDevices colorRampPalette
#' @importFrom stats quantile
#' @noRd

.makeLegend <- function(scaleValues, colours = NULL, title = NULL,
                        prevX = unit(0, "points"), theme){

  legends <- list()
  allLabels <- scaleValues
  if(is.null(colours)) colours <- theme$colours

  if(length(allLabels) > 10){
    testItems <- sample(allLabels, 10)
  } else {
    testItems <- allLabels
  }
  if(any(testItems %in% grDevices::colors() | any(grepl(pattern = "\\#(.{6,8})", x = testItems)))){
    return(legends)
  }

  if(is.null(allLabels)){
    return(legends)
  }

  # the ramp must have exactly one entry per scale value: .makeScale has already
  # decided how many colour steps there are, and the raster is coloured from that
  # same vector. Deriving the count from anything else here would let the legend
  # claim a colour the map never draws.
  thebins <- length(allLabels)

  # determine the tick values and labels
  if(thebins > theme$legend.bins){
    tickPositions <- quantile(1:thebins, probs = seq(0, 1, length.out = theme$legend.bins + 1), type = 1, names = FALSE)
  } else {
    tickPositions <- 1:thebins
  }
  legendLabels <- allLabels[tickPositions]

  if(!theme$legend.ascending){
    tickPositions <- rev(tickPositions)
  }

  legend_values <- textGrob(label = legendLabels,
                            name = "legend_values",
                            gp = gpar(col = NA))

  # colour ramp legend
  allColours <- .resolve_colours(colours, thebins)
  if(!is.null(names(allColours))) allColours <- unname(allColours)

  legend_obj <- rasterGrob(x = unit(0, "npc") + unit(5, "points") + prevX,
                           width = unit(10, "points"),
                           height = unit(1, "npc"),
                           just = c("left"),
                           name = "legend_items",
                           image = rev(allColours),
                           interpolate = FALSE)

  if(theme$legend.box.plot){
    legend_obj <- gList(
      legend_obj,
      rectGrob(x = unit(0, "npc") + unit(5, "points") + prevX,
               width = unit(10, "points"),
               just = c("left"),
               name = "legend_box",
               gp = gpar(col = theme$legend.box.colour,
                         fill = NA,
                         lty = theme$legend.box.linetype,
                         lwd = theme$legend.box.linewidth)))
  }

  legend_labels <- NULL
  if(theme$legend.label.plot){
    thePositions <- tickPositions / max(tickPositions)
    thePositions <- (tickPositions - 1) / max(tickPositions) + thePositions[1] / 2

    if(is.numeric(legendLabels)){
      # rounded first, so a value next to 0 reads 0.00 rather than 2.6e-06
      legendLabels <- format(round(legendLabels, theme$legend.digits + 1),
                             digits = theme$legend.digits + 1,
                             scientific = FALSE)
    }

    legend_labels <- textGrob(label = legendLabels,
                              x = unit(0, "npc") + unit(20, "points") + prevX,
                              y = unit(thePositions, "npc"),
                              name = "legend_labels",
                              just = c("left", "centre"),
                              gp = gpar(fontsize = theme$legend.label.fontsize,
                                        col = theme$legend.label.colour))
  }

  # a title is only meaningful when several scales share the panel, but drawing
  # it whenever asked keeps the rule simple
  legend_title <- NULL
  if(!is.null(title)){
    legend_title <- textGrob(label = title,
                             x = unit(0, "npc") + unit(5, "points") + prevX,
                             y = unit(1, "npc") + unit(7, "points"),
                             name = "legend_title",
                             just = c("left", "bottom"),
                             gp = gpar(fontsize = theme$legend.label.fontsize,
                                       col = theme$legend.label.colour))
  }

  out <- gTree(children = gList(legend_values, legend_obj, legend_labels,
                                legend_title))
  legends <- c(legends, list(fillcol = out))

  return(legends)
}


#' Make the raster grob of a plot
#'
#' Grid-only: creates a \code{rasterGrob} from plot values.
#' @param plotValues numeric vector of cell values.
#' @param scaleValues sorted unique values for the colour scale.
#' @param rows number of rows in the grid.
#' @param cols number of columns in the grid.
#' @param layout the layout object from \code{.makeLayout}.
#' @param colours the layer's colour specification, or NULL for the theme's.
#' @param categories the layer's category table, used to match a named colour
#'   vector to labels.
#' @param theme the completed theme list.
#' @importFrom grDevices colorRampPalette
#' @importFrom grid gpar unit rasterGrob gList
#' @importFrom checkmate testCharacter
#' @noRd

.makeGrob <- function(plotValues, scaleValues, rows, cols, layout,
                      colours = NULL, categories = NULL, theme){

  if(is.null(theme)) theme <- .themeDefaults
  if(is.null(colours)) colours <- theme$colours

  if(!theme$box.plot) return(NULL)

  if(testCharacter(x = plotValues, pattern = "\\#(.{6,8})")){
    theColours <- as.vector(plotValues)

  } else if(!is.null(names(colours))){
    # named vector: an explicit category -> colour map. Matched by the category
    # label where the layer has a category table, else by the raw cell value,
    # so a layer without categories can still be mapped by its gids.
    theColours <- rep(NA_character_, length(plotValues))
    if(!is.null(categories) && !is.null(categories$gid)){
      lab <- categories$val
      if(is.null(lab)) lab <- as.character(categories$gid)
      known <- names(colours)
      unmatched <- setdiff(lab, known)
      if(length(unmatched) > 0){
        stop("no colour given for categor",
             if(length(unmatched) == 1) "y " else "ies ",
             paste0("'", unmatched, "'", collapse = ", "),
             ". A named colour vector must cover every category present, ",
             "otherwise the unnamed ones would be drawn on a different scale.",
             call. = FALSE)
      }
      for(k in seq_along(categories$gid)){
        theColours[!is.na(plotValues) & plotValues == categories$gid[k]] <-
          colours[[lab[k]]]
      }
    } else {
      for(k in seq_along(colours)){
        theColours[!is.na(plotValues) &
                     as.character(plotValues) == names(colours)[k]] <- colours[[k]]
      }
    }

  } else if(length(scaleValues) <= 1){
    # zero or one unique value: uniform colour for non-NA cells
    allColours <- .resolve_colours(colours, max(length(scaleValues), 1))
    theColours <- rep(NA_character_, length(plotValues))
    theColours[!is.na(plotValues)] <- allColours[1]

  } else {
    allColours <- .resolve_colours(colours, length(scaleValues))

    # clamp to the scale's ends: with an explicit limits/shared scale a cell
    # may lie outside the breaks, and cut() would return NA for it -- a hole
    # in the map that reads as missing data rather than as out-of-range.
    clamped <- plotValues
    if(is.numeric(clamped)){
      lo <- scaleValues[1]
      hi <- scaleValues[length(scaleValues)]
      clamped[!is.na(clamped) & clamped < lo] <- lo
      clamped[!is.na(clamped) & clamped > hi] <- hi
    }

    scaleBreaks <- c(scaleValues[1] - 1, scaleValues)
    valCuts <- cut(clamped, breaks = scaleBreaks, include.lowest = TRUE)
    theColours <- allColours[valCuts]
  }

  if(!is.null(theme$colours.missing)){
    theColours[is.na(theColours)] <- theme$colours.missing
  }

  winX <- layout$window$xmin
  winY <- layout$window$ymin
  winW <- layout$window$xmax - layout$window$xmin
  winH <- layout$window$ymax - layout$window$ymin

  out <- rasterGrob(x = unit(winX - layout$margin$x, "native"),
                    y = unit(winY - layout$margin$y, "native"),
                    width = unit(winW + 2 * layout$margin$x, "native") - unit(theme$box.linewidth, "points"),
                    height = unit(winH + 2 * layout$margin$y, "native") - unit(theme$box.linewidth, "points"),
                    hjust = 0,
                    vjust = 0,
                    image = matrix(data = theColours, nrow = rows, ncol = cols, byrow = TRUE),
                    name = "theRaster",
                    interpolate = FALSE)

  gList(out)
}


#' Derive the colour scale for a layer
#'
#' Decides, once, the values the colour ramp is built over. Both the raster and
#' the legend index this same vector, so they cannot disagree about what a
#' colour means.
#'
#' Two regimes:
#' \itemize{
#'   \item \strong{Categorical} -- few distinct values, or non-numeric. Every
#'     distinct value gets its own colour, which is the historical behaviour and
#'     the right one for class maps.
#'   \item \strong{Continuous} -- more distinct values than \code{bins}. The
#'     range is cut into \code{bins} equal steps. Without this a layer with
#'     50,000 distinct elevations builds a 50,000-entry ramp per panel, which is
#'     both slow and gives the legend one label per raw value.
#' }
#'
#' \code{limits} overrides the range taken from the data. Passing the same
#' \code{limits} to several panels is what makes them comparable: identical
#' values then map to identical colours regardless of each panel's own spread.
#' Values outside \code{limits} are clamped to the end colours rather than
#' dropped, so an out-of-range cell stays visible instead of turning into a
#' hole.
#'
#' @param values numeric (or other) vector of cell values, NAs already removed.
#' @param bins [`numeric(1)`][numeric]\cr maximum number of colour steps before
#'   the scale switches from categorical to binned-continuous.
#' @param limits [`numeric(2)`][numeric]\cr \code{c(min, max)} to force, or
#'   \code{NULL} to take the range from \code{values}.
#' @return numeric vector of scale values (ascending), or the sorted unique
#'   values when the layer is categorical.
#' @noRd

.makeScale <- function(values, bins = 256, limits = NULL){

  if(length(values) == 0) return(numeric(0))

  if(!is.numeric(values)){
    return(sort(unique(values)))
  }

  uniq <- as.numeric(sortUniqueCpp(values))

  # categorical: keep one colour per distinct value unless limits force a
  # continuous reading of the layer
  if(is.null(limits) && length(uniq) <= bins) return(uniq)

  if(is.null(limits)){
    lo <- uniq[1]
    hi <- uniq[length(uniq)]
  } else {
    lo <- limits[1]
    hi <- limits[2]
  }

  # a degenerate range has no interior to cut
  if(!is.finite(lo) || !is.finite(hi) || hi <= lo) return(lo)

  seq(from = lo, to = hi, length.out = bins)
}


#' Compute a hillshade multiplier from an elevation layer
#'
#' Standard analytical hillshade: the cosine of the angle between the surface
#' normal and a directional light source, following the same formulation ArcGIS
#' and GDAL use. Slope and aspect come from \code{slopeAspectCpp} (Horn's
#' method), so the estimator is shared with \code{make_radiation} and cannot
#' drift from it.
#'
#' The returned value is an illumination factor in \[0, 1\], not a colour: 1 is
#' fully lit, 0 is fully shadowed. Callers multiply it into the base colours.
#' This is local (self-shading) hillshade only -- it does not cast shadows from
#' distant terrain, which would need a ray march per cell.
#'
#' Elevation is assumed to be in the same linear unit as the cell size (metres,
#' per the absolute-units convention). \code{exaggeration} scales the elevation
#' before the gradient is taken, which is how a subtle relief is made legible
#' without altering the data.
#'
#' @param obj a mosaik object.
#' @param layer [`character(1)`][character]\cr name of the elevation layer.
#' @param azimuth [`numeric`][numeric]\cr direction of the light source in degrees
#'   clockwise from north. Default 315 (north-west), the cartographic convention
#'   -- light from the upper left avoids the relief-inversion illusion that
#'   lighting from below triggers in most viewers. Several azimuths average into
#'   a multi-directional hillshade, which keeps slopes facing away from a single
#'   light from going flat and featureless.
#' @param altitude [`numeric(1)`][numeric]\cr height of the light source in
#'   degrees above the horizon. Default 45.
#' @param exaggeration [`numeric(1)`][numeric]\cr vertical exaggeration applied to
#'   the elevation before the gradient. Default 1 (no exaggeration).
#' @return numeric vector of illumination factors in \[0, 1\], row-major, one per
#'   cell; NA where the elevation is NA.
#' @importFrom checkmate assertNumber assertChoice
#' @noRd

.hillshade <- function(obj, layer, azimuth = 315, altitude = 45,
                       exaggeration = 1){

  assertNumeric(x = azimuth, lower = 0, upper = 360, min.len = 1,
                any.missing = FALSE)
  assertNumber(x = altitude, lower = 0, upper = 90)
  assertNumber(x = exaggeration, lower = 0, finite = TRUE)

  elev <- as.numeric(msk_pull(obj, layer))
  # slope/aspect/TWI rules assume elevation in metres; a raw drw_texture /
  # drw_noise layer lies in [0,1] and yields degenerate slope with no error.
  # A genuinely flat sub-metre DEM also lies in [0,1] and cannot be told
  # apart by range alone, hence the escape hatch in the message.
  .rng <- range(elev, na.rm = TRUE)
  if (.rng[1] >= 0 && .rng[2] <= 1) {
    warning("surface layer '", layer, "' has values within [0, 1] (range ",
            round(.rng[1], 3), "-", round(.rng[2], 3), "), which looks like ",
            "unscaled pattern output. Slope and TWI assume elevation in ",
            "metres -- mdf_scale() it to a real relief (e.g. range = c(0, 800)) ",
            "if so. Ignore this warning if the [0, 1] scale is deliberate.",
            call. = FALSE)
  }

  nrows <- obj@dims[2]
  ncols <- obj@dims[1]
  res <- msk_res(obj)

  # exaggerate by scaling elevation, not by scaling the resulting slope: the
  # two differ once slope is non-linear in the gradient (atan), and scaling the
  # input is what "vertical exaggeration" means.
  sa <- slopeAspectCpp(elev = elev * exaggeration, nrow = nrows, ncol = ncols,
                       xres = res[1], yres = res[2])
  slope  <- sa$slope      # radians from horizontal
  aspect <- sa$aspect     # radians clockwise from north, -1 = flat

  zenith <- (90 - altitude) * pi / 180

  # flat cells have no aspect; the aspect term vanishes there anyway (slope 0),
  # but aspect = -1 would otherwise feed a meaningless angle into the cosine.
  aspect[!is.na(aspect) & aspect < 0] <- 0

  # one light per azimuth, averaged. The slope/aspect pass above is shared, so
  # additional lights cost only the cosine.
  shade <- numeric(length(slope))
  for(a in azimuth){
    az <- a * pi / 180
    s <- cos(zenith) * cos(slope) +
      sin(zenith) * sin(slope) * cos(az - aspect)
    # self-shadowed faces clamp to 0 before averaging, not after: a face turned
    # away from one light contributes nothing to that light, rather than
    # subtracting from the others.
    s[!is.na(s) & s < 0] <- 0
    shade <- shade + s
  }
  shade <- shade / length(azimuth)

  shade[!is.na(shade) & shade > 1] <- 1

  shade
}


#' Multiply base colours by a hillshade factor
#'
#' Composites the illumination onto already-resolved colours in RGB space.
#' \code{intensity} mixes between the unshaded colour (0) and the fully shaded
#' one (1), so the relief can be dialled back without recomputing it.
#'
#' @param colours character vector of hex colours (may contain NA).
#' @param shade numeric vector of illumination factors in \[0, 1\].
#' @param intensity [`numeric(1)`][numeric]\cr strength of the effect in \[0, 1\].
#' @return character vector of hex colours, NA preserved.
#' @importFrom grDevices col2rgb rgb
#' @noRd

.apply_hillshade <- function(colours, shade, intensity = 1){

  keep <- !is.na(colours) & !is.na(shade)
  if(!any(keep)) return(colours)

  # blend the multiplier toward 1 so intensity < 1 lightens the whole effect
  f <- 1 - intensity + intensity * shade[keep]

  rgbMat <- col2rgb(colours[keep], alpha = TRUE)
  out <- colours
  out[keep] <- rgb(red = pmin(255, rgbMat["red", ] * f),
                   green = pmin(255, rgbMat["green", ] * f),
                   blue = pmin(255, rgbMat["blue", ] * f),
                   alpha = rgbMat["alpha", ],
                   maxColorValue = 255)
  out
}


#' Prepare one panel for plotting
#'
#' Builds the raster and the legends for a panel, compositing its layers bottom
#' to top. Each layer gets its own colour scale and its own legend, so a
#' variable drawn over terrain is read against its own numbers rather than the
#' terrain's.
#'
#' @param obj a mosaik object.
#' @param specs list of \code{mskLayer} specifications for this panel.
#' @param sharedLimits named list of pooled ranges, by layer name.
#' @param window named list with xmin, xmax, ymin, ymax.
#' @param theme the completed theme list.
#' @return A list with elements \code{theme}, \code{grob}, \code{legend},
#'   \code{layout}.
#' @importFrom grid unit convertX
#' @noRd

.makePlot <- function(obj, specs, sharedLimits = list(), window, theme){

  if(is.null(theme)) theme <- .themeDefaults

  out <- list(theme = theme, grob = NULL, legend = NULL, layout = NULL)

  rows <- obj@dims[2]
  cols <- obj@dims[1]

  baseMat <- NULL
  theGrob <- NULL
  theLegend <- list()
  prevX <- unit(0, "points")

  for(s in specs){

    plotValues <- msk_pull(obj, s$layer)
    nonNA <- plotValues[!is.na(plotValues)]

    # limits: the layer's own beats the shared pool beats the theme
    theLimits <- s$limits
    if(is.null(theLimits)) theLimits <- sharedLimits[[s$layer]]
    if(is.null(theLimits)) theLimits <- theme$scale.range

    theBins <- s$bins
    if(is.null(theBins)) theBins <- theme$scale.bins
    if(is.null(theBins)) theBins <- 256

    scaleValues <- .makeScale(values = nonNA, bins = theBins,
                              limits = theLimits)

    theColours <- s$colours
    # a category table may carry a colour per class (as landscape's land cover
    # does); it is used unless the call gives colours
    cats <- obj@categories[[s$layer]]
    if(is.null(theColours) && !is.null(cats$colour)){
      theColours <- stats::setNames(cats$colour,
                                    if(is.null(cats$val)) cats$gid else cats$val)
    }
    if(is.null(theColours)) theColours <- theme$colours

    # layout needs the legends to size the plotting area, but the legend needs
    # nothing from the layout, so legends are built first
    if(theme$legend.plot && isTRUE(s$legend)){
      theTitle <- s$title
      # only title the strips when several share the panel: a lone legend needs
      # no heading, the panel title already says what is drawn
      if(is.null(theTitle) && length(specs) > 1) theTitle <- s$layer
      # a named colour vector may hold classes the map does not show; the
      # legend lists the values shown, so it takes their colours by name
      legendColours <- theColours
      if(!is.null(names(theColours))){
        key <- as.character(scaleValues)
        if(!is.null(cats$gid) && !is.null(cats$val)){
          key <- cats$val[match(scaleValues, cats$gid)]
        }
        if(!anyNA(theColours[key])) legendColours <- unname(theColours[key])
      }
      thisLegend <- .makeLegend(scaleValues = scaleValues,
                                colours = legendColours,
                                title = theTitle,
                                prevX = prevX,
                                theme = theme)
      if(length(thisLegend) > 0){
        theLegend <- c(theLegend, thisLegend)
        # offset the next strip past this one's widest label
        labels <- thisLegend[[1]]$children$legend_labels$label
        if(!is.null(labels)){
          maxLbl <- labels[which.max(nchar(labels))]
          prevX <- prevX + unit(1, "strwidth", maxLbl) + unit(25, "points")
        } else {
          prevX <- prevX + unit(35, "points")
        }
      }
    }

    theLayout <- .makeLayout(obj = obj, legend = theLegend,
                             window = window, theme = theme)

    thisGrob <- .makeGrob(plotValues = plotValues,
                          scaleValues = scaleValues,
                          rows = rows, cols = cols,
                          layout = theLayout,
                          colours = theColours,
                          categories = obj@categories[[s$layer]],
                          theme = theme)
    if(is.null(thisGrob)) next

    # hillshade multiplies this layer's own colours before it is composited, so
    # a layer drawn on top of a shaded one is not darkened by that shading
    if(!is.null(s$hillshade)){
      if(!s$hillshade$layer %in% names(obj@layers)){
        warning("hillshade layer '", s$hillshade$layer, "' does not exist; ",
                "skipping the hillshade.", call. = FALSE)
      } else {
        shadeArgs <- s$hillshade
        intensity <- if(is.null(shadeArgs$intensity)) 1 else shadeArgs$intensity
        shadeArgs$intensity <- NULL
        shadeArgs$layer <- NULL

        shade <- do.call(.hillshade,
                         c(list(obj = obj, layer = s$hillshade$layer),
                           shadeArgs))

        # a grid raster object flattens by ROW, so `as.vector()` returns the
        # cells in the same order as the shade. The write-back has to be a
        # raster object again: assigning a plain matrix here drew every shaded
        # panel transposed, because grid reads a bare matrix by column.
        shaded <- .apply_hillshade(colours = as.vector(thisGrob[[1]]$raster),
                                   shade = shade, intensity = intensity)
        thisGrob[[1]]$raster <- grDevices::as.raster(matrix(shaded, nrow = rows,
                                                 ncol = cols, byrow = TRUE))
      }
    }

    # composite onto what is already there: NA is transparent, so the layer
    # below shows through wherever this one has no value
    if(is.null(baseMat)){
      baseMat <- thisGrob[[1]]$raster
      theGrob <- thisGrob
    } else {
      top <- thisGrob[[1]]$raster
      drawn <- !is.na(top) & top != theme$colours.missing
      baseMat[drawn] <- top[drawn]
      theGrob[[1]]$raster <- baseMat
    }
  }

  out$grob <- theGrob
  out$legend <- theLegend
  out$layout <- .makeLayout(obj = obj, legend = theLegend,
                            window = window, theme = theme)

  return(out)
}

