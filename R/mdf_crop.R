#' Crop a mosaik to a smaller extent
#'
#' Subset a mosaik to the cells that fall within a new extent. The output grid
#' aligns to the input resolution — the requested extent is snapped outward to
#' cell boundaries.
#' @param obj [mosaik]\cr the mosaik to crop.
#' @param extent [numeric(4)][numeric]\cr the target extent as
#'   \code{c(xmin, xmax, ymin, ymax)}. Values outside the input extent are
#'   silently clamped.
#' @return A mosaik with reduced extent and dimensions. Layer values are
#'   subsetted; categories, patches, and global are preserved (patches/global
#'   may reference values no longer present).
#' @family utilities
#' @importFrom checkmate assertClass assertNumeric
#' @export

mdf_crop <- function(obj = NULL,
                     extent) {

  if (.is_recipe(obj)) return(.record_step(obj, match.call()))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumeric(x = extent, len = 4, any.missing = FALSE)
  if (extent[2] <= extent[1]) stop("extent xmax must be > xmin.")
  if (extent[4] <= extent[3]) stop("extent ymax must be > ymin.")

  old_ext  <- obj@extent
  dims     <- obj@dims   # c(ncols, nrows)
  res      <- msk_res(obj)

  # clamp requested extent to input extent
  xmin <- max(extent[1], old_ext[1])
  xmax <- min(extent[2], old_ext[2])
  ymin <- max(extent[3], old_ext[3])
  ymax <- min(extent[4], old_ext[4])

  if (xmin >= xmax || ymin >= ymax) {
    stop("Requested extent does not overlap with the mosaik extent.")
  }

  # snap to cell boundaries (outward)
  col_start <- floor((xmin - old_ext[1]) / res[1]) + 1L
  col_end   <- ceiling((xmax - old_ext[1]) / res[1])
  row_start <- floor((old_ext[4] - ymax) / res[2]) + 1L
  row_end   <- ceiling((old_ext[4] - ymin) / res[2])

  # clamp to valid range
  col_start <- max(1L, as.integer(col_start))
  col_end   <- min(dims[1], as.integer(col_end))
  row_start <- max(1L, as.integer(row_start))
  row_end   <- min(dims[2], as.integer(row_end))

  new_ncols <- col_end - col_start + 1L
  new_nrows <- row_end - row_start + 1L

  # compute snapped extent from cell boundaries
  new_ext <- c(
    old_ext[1] + (col_start - 1L) * res[1],
    old_ext[1] + col_end * res[1],
    old_ext[4] - row_end * res[2],
    old_ext[4] - (row_start - 1L) * res[2]
  )

  # subset layers ----
  new_layers <- list()
  for (nm in names(obj@layers)) {
    vals <- msk_pull(obj, nm)
    # vals is stored row-major: row 1 = cells 1..ncols, row 2 = ncols+1..2*ncols
    idx <- unlist(lapply(row_start:row_end, function(r) {
      ((r - 1L) * dims[1]) + col_start:col_end
    }))
    new_layers[[nm]] <- vals[idx]
  }

  # build output ----
  prov <- msk_prov("mdf_crop", list(extent = extent))

  new_mosaik(
    extent     = new_ext,
    dims       = c(new_ncols, new_nrows),
    layers     = new_layers,
    categories = obj@categories,
    patches    = obj@patches,
    global     = obj@global,
    crs        = obj@crs,
    provenance = c(obj@provenance, list(prov))
  )
}
