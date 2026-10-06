#' Summarise a layer within or around each zone
#'
#' Group cells by the zone they belong to and replace every cell of a zone with
#' a single summary of another layer, computed either over the zone's own cells
#' or over the cells adjacent to it. This is the raster counterpart of
#' \code{summarise(.by = )} in dplyr: the zones are the groups and \code{layer}
#' supplies the values. Where dplyr returns one row per group, a layer needs a
#' value in every cell, so the summary is written onto every cell of its zone.
#'
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param by [`character(1)`][character]\cr the zones: the name of a layer of
#'   labels in \code{obj} (e.g. from \code{\link{mdf_componentise}}). Cells that
#'   are \code{0} or \code{NA} are outside every zone. To use zones held in
#'   another mosaik, bring them in with \code{\link{msk_add}} first.
#' @param fun [`character(1)`][character] or [`function`][function]\cr how to
#'   summarise the values of a zone. Either a function taking a numeric vector
#'   and returning a single value, or one of the shorthands \code{"max"}, \code{"min"}, \code{"sum"},
#'   \code{"mean"}, \code{"median"}, \code{"any"} (any non-zero value present),
#'   \code{"all"}, \code{"n"} (number of cells), \code{"n_distinct"} (number of
#'   distinct non-zero values) or \code{"unique"} (the one non-zero value
#'   present, \code{NA} if there are several, e.g. the patch around a hole).
#' @param neighbours [`logical(1)`][logical]\cr what the summary is taken over. If
#'   \code{FALSE} (default), the zone's own cells. If \code{TRUE}, the cells
#'   adjacent to the zone but not part of it, which answers questions of the
#'   form "what does this zone touch".
#' @param connectivity [`integerish(1)`][integer]\cr which cells count as adjacent
#'   when \code{neighbours = TRUE}: \code{8} (default, queen) or \code{4} (rook).
#'   Ignored otherwise.
#' @param background [`numeric(1)`][numeric]\cr the value written to cells that
#'   belong to no zone. Defaults to \code{NA}.
#' @param layer [`character(1)`][character]\cr the layer holding the values to
#'   summarise. Defaults to the first layer. Use the same layer as \code{by} to
#'   summarise the zones by their own labels.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimensions as \code{obj}, where every cell of a
#'   zone carries that zone's summary value.
#' @details
#'   \code{mdf_summarise} answers questions that need a whole connected object to be
#'   judged at once and the verdict handed back to each of its cells. No
#'   moving-window operator can do this, because the object's extent is not
#'   known in advance. Two shapes of question, selected by \code{neighbours}:
#'
#'   \describe{
#'     \item{\code{neighbours = FALSE}}{summarise the zone's own cells. "Does
#'       this patch contain any core cell?" is \code{fun = "any"} over a core
#'       mask; "how far does this run reach?" is \code{fun = "max"} over a
#'       distance layer.}
#'     \item{\code{neighbours = TRUE}}{summarise what surrounds the zone. "How
#'       many distinct patches does this connector touch?" is
#'       \code{fun = "n_distinct"} over a patch-ID layer. Paired with
#'       \code{fun = "n"} for the zone's own size, that is what separates a
#'       bridge from a loop in MSPA: a run touching two or more cores is a
#'       bridge only if it is more than one cell long, since a single cell
#'       wedged diagonally between two cores does not cross between them.}
#'   }
#'
#'   Zones need not be spatially connected: any label layer works, so classes,
#'   administrative units or strata are equally valid groups. When zones are
#'   connected components, pair it with \code{\link{mdf_componentise}}.
#' @examples
#' # the size of each forest patch, written onto every cell of the patch
#' m <- landscape |>
#'   mdf_filter(cover == 47, add = "forest") |>
#'   mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
#'   mdf_summarise(by = "patch", fun = "n", layer = "patch", add = "size")
#' msk_vis(m, .layer("forest"), .layer("size"))
#'
#' # how many distinct forest patches each open area touches
#' m <- m |>
#'   mdf_filter(forest == 0, add = "open") |>
#'   mdf_componentise(connectivity = 4L, layer = "open", add = "area") |>
#'   mdf_summarise(by = "area", fun = "n_distinct", neighbours = TRUE,
#'                 layer = "patch", add = "touching")
#' msk_vis(m, .layer("forest"), .layer("open"), .layer("touching"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertCharacter assertLogical assertChoice
#'   assertNumber
#' @export

mdf_summarise <- function(obj = NULL,
                      by = NULL,
                      fun = "sum",
                      neighbours = FALSE,
                      connectivity = 8L,
                      background = NA,
                      layer = NULL,
                      add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertLogical(x = neighbours, len = 1, any.missing = FALSE)
  assertChoice(x = connectivity, choices = c(4L, 8L))
  assertNumber(x = background, na.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)
  if(is.null(by)) stop("'by' is required: the layer holding the zones.")

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  zones <- .pull_layer(obj, by, "by")
  dims <- obj@dims
  ncols <- dims[1]; nrows <- dims[2]

  # resolve the summary function ----
  f <- .summary_fun(fun, "mdf_summarise")

  # body ----
  # cells belonging to a zone; 0/NA is "no zone" and is never summarised
  inZone <- !is.na(zones) & zones != 0
  temp <- rep(as.numeric(background), length(vals))

  if(any(inZone)){
    zoneIds <- zones[inZone]
    cellIdx <- which(inZone)

    if(!neighbours){
      # summarise each zone over its own cells
      groups <- split(cellIdx, zoneIds)
      for(zi in names(groups)){
        temp[groups[[zi]]] <- f(vals[groups[[zi]]])
      }

    } else {
      # summarise each zone over the cells adjacent to it but outside it. The
      # neighbour set is gathered by shifting the whole grid, which keeps the
      # cost proportional to the grid rather than to the number of zones.
      zMat <- matrix(zones, nrow = nrows, ncol = ncols, byrow = TRUE)
      vMat <- matrix(vals,  nrow = nrows, ncol = ncols, byrow = TRUE)

      offsets <- if(connectivity == 8L){
        list(c(-1,-1), c(-1,0), c(-1,1), c(0,-1), c(0,1), c(1,-1), c(1,0), c(1,1))
      } else {
        list(c(-1,0), c(0,-1), c(0,1), c(1,0))
      }

      # for every shift, pair each in-zone cell with the value of the cell it
      # looks at, keeping only pairs that leave the zone
      pairsZone <- integer(0)
      pairsVal  <- numeric(0)
      for(off in offsets){
        dr <- off[1]; dc <- off[2]
        rSrc <- max(1, 1 - dr):min(nrows, nrows - dr)
        cSrc <- max(1, 1 - dc):min(ncols, ncols - dc)
        if(length(rSrc) == 0 || length(cSrc) == 0) next

        zSrc <- zMat[rSrc, cSrc, drop = FALSE]
        zNbr <- zMat[rSrc + dr, cSrc + dc, drop = FALSE]
        vNbr <- vMat[rSrc + dr, cSrc + dc, drop = FALSE]

        keep <- !is.na(zSrc) & zSrc != 0 &
                (is.na(zNbr) | zNbr != zSrc)   # the neighbour is outside the zone
        if(!any(keep)) next

        pairsZone <- c(pairsZone, zSrc[keep])
        pairsVal  <- c(pairsVal,  vNbr[keep])
      }

      if(length(pairsZone) > 0){
        summaries <- vapply(split(pairsVal, pairsZone), f, numeric(1))
        # map each zone's summary back onto its cells
        hit <- match(as.character(zoneIds), names(summaries))
        temp[cellIdx] <- unname(summaries[hit])
      }
    }
  }

  # build output ----
  # a zonal summary is a new kind of value, not the input's
  .update_mosaik(obj, values = temp, keep = FALSE, step = step)
}
