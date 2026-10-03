#' Utility functions that act on the tabular information of a mosaik
#'
#' These utility functions are meant to immitate the tidy dplyr logic on the
#' tables of a mosaik. \code{msk_select} keeps only the specified layers,
#' \code{msk_remove} drops them, \code{msk_add} brings layers in from a second
#' mosaik on the same grid, or from vectors of cell values. \code{msk_pull}
#' extracts cell values from a layer as a vector (like \code{dplyr::pull}).
#'
#' @param obj [`mosaik`][mosaik]\cr the mosaik object.
#' @param ... layer names (unquoted or character strings). For
#'   \code{msk_select} and \code{msk_remove} these name layers of \code{obj},
#'   for \code{msk_add} they name layers of \code{from}. \code{msk_add} also
#'   takes \code{name = values} pairs, which add a vector of cell values as a
#'   new layer of that name.
#' @name utils
NULL

# --- msk_select ---------------------------------------------------------------

#' @rdname utils
#' @return \code{msk_select}: A mosaik containing only the selected layers.
#' @examples
#' # keep one layer, dropping everything else
#' msk_select(landscape, cover)
#'
#' # layer names may also be given as strings
#' msk_select(landscape, "cover", "canopy")
#' @importFrom checkmate assertClass
#' @export

msk_select <- function(obj, ...){

  step <- .step()
  assertClass(x = obj, classes = "mosaik")

  vars <- as.character(match.call(expand.dots = FALSE)$...)

  layerNames <- names(obj@layers)
  keep <- vars[vars %in% layerNames]
  if(length(keep) == 0){
    stop("none of the requested layers found: ", paste(vars, collapse = ", "))
  }

  .update_mosaik(obj,
                 layers = obj@layers[keep],
                 categories = obj@categories[intersect(keep, names(obj@categories))],
                 patches = .unlink_patches(obj@patches[intersect(keep, names(obj@patches))],
                                           keep),
                 global = obj@global[intersect(keep, names(obj@global))],
                 step = step)
}

# --- msk_remove ---------------------------------------------------------------

#' @rdname utils
#' @return \code{msk_remove}: A mosaik without the named layers.
#' @examples
#' # drop a layer, keeping the rest
#' msk_remove(landscape, canopy)
#'
#' # removing every layer is refused
#' try(msk_remove(landscape, cover, canopy))
#' @importFrom checkmate assertClass
#' @export

msk_remove <- function(obj, ...){

  step <- .step()
  assertClass(x = obj, classes = "mosaik")

  vars <- as.character(match.call(expand.dots = FALSE)$...)

  layerNames <- names(obj@layers)
  drop <- vars[vars %in% layerNames]
  keep <- setdiff(layerNames, drop)
  if(length(keep) == 0){
    stop("removing ", paste(drop, collapse = ", "), " would leave no layer.")
  }

  .update_mosaik(obj,
                 layers = obj@layers[keep],
                 categories = obj@categories[intersect(keep, names(obj@categories))],
                 patches = .unlink_patches(obj@patches[intersect(keep, names(obj@patches))],
                                           keep),
                 global = obj@global[intersect(keep, names(obj@global))],
                 step = step)
}

# --- msk_add ------------------------------------------------------------------

#' @rdname utils
#' @param from [`mosaik`][mosaik]\cr the mosaik to take layers from. Must sit on
#'   the same grid as \code{obj}, i.e. carry the same \code{extent},
#'   \code{dims} and \code{crs}. Not needed when only vectors of values are
#'   added.
#' @param rename [`character(.)`][character]\cr optional names under which the
#'   layers of \code{from} should be stored in \code{obj}, in the order they
#'   are given in \code{...}. Needed when a name is already taken.
#' @return \code{msk_add}: A mosaik with the layers added.
#' @details \code{msk_add} is the inverse of \code{msk_select} and the way to
#'   combine two mosaiks: operators such as \code{\link{mdf_blend}} work within
#'   a single object, so a layer computed elsewhere is brought in first and
#'   combined afterwards. The grid is checked here, which is why no operator
#'   needs to check it again. Whatever was measured on an added layer (its
#'   class, patch and landscape-level results) travels with it. Its patches
#'   stay measurable only if the layer holding their numbers (see
#'   \code{\link{mdf_componentise}}) is added as well.
#' @examples
#' # a second mosaik on the same grid as 'landscape'
#' other <- mosaik(extent = c(0, 60, 0, 56), res = 1,
#'                 vals = list(elevation = runif(60 * 56, 0, 800)))
#'
#' # bring one layer over, then use it alongside the layers already there
#' combined <- msk_add(landscape, other, elevation)
#' msk_names(combined)
#'
#' # a name already taken must be renamed
#' twin <- msk_select(landscape, cover)
#' try(msk_add(landscape, twin, cover))
#' msk_names(msk_add(landscape, twin, cover, rename = "cover2"))
#'
#' # a vector of cell values becomes a layer of its own
#' seed <- as.numeric(seq_len(msk_ncells(landscape)) == 100)
#' msk_names(msk_add(landscape, seed = seed))
#' @importFrom checkmate assertClass assertCharacter
#' @export

msk_add <- function(obj, from = NULL, ..., rename = NULL){

  step <- .step()
  assertClass(x = obj, classes = "mosaik")
  assertClass(x = from, classes = "mosaik", null.ok = TRUE)
  assertCharacter(x = rename, null.ok = TRUE, any.missing = FALSE)

  # unnamed arguments name layers of 'from', named ones are vectors of values
  dots <- match.call(expand.dots = FALSE)$...
  nms <- names(dots)
  if(is.null(nms)) nms <- rep("", length(dots))
  vars <- as.character(dots[!nzchar(nms)])
  values <- list()
  for(i in which(nzchar(nms))) values[[nms[i]]] <- ...elt(i)

  new_layers <- obj@layers
  new_categories <- obj@categories
  new_patches <- obj@patches
  new_global <- obj@global
  provenance <- obj@provenance

  if(!is.null(from)){

    # grid conformity ----
    if(!isTRUE(all.equal(obj@extent, from@extent))){
      stop("'from' does not sit on the same extent as 'obj'.", call. = FALSE)
    }
    if(!identical(obj@dims, from@dims)){
      stop("'from' does not have the same dimensions as 'obj'.", call. = FALSE)
    }
    if(!identical(obj@crs, from@crs)){
      stop("'from' does not have the same crs as 'obj'.", call. = FALSE)
    }

    # resolve which layers to take ----
    if(length(vars) == 0 && length(values) == 0) vars <- names(from@layers)

    missing <- setdiff(vars, names(from@layers))
    if(length(missing) > 0){
      stop("layer(s) not found in 'from': ", paste(missing, collapse = ", "),
           call. = FALSE)
    }

    # resolve target names ----
    if(is.null(rename)){
      target <- vars
    } else {
      if(length(rename) != length(vars)){
        stop("'rename' must name as many layers as are added.", call. = FALSE)
      }
      target <- rename
    }

    for(i in seq_along(vars)){
      new_layers[[target[i]]] <- from@layers[[vars[i]]]
      cats <- from@categories[[vars[i]]]
      if(!is.null(cats)) new_categories[[target[i]]] <- cats
      # the patch record keeps its link to the layer holding the patch numbers
      # only if that layer comes along too, under its new name
      rec <- from@patches[[vars[i]]]
      if(!is.null(rec)){
        if(!is.null(rec$ids)){
          j <- match(rec$ids, vars)
          if(is.na(j)) rec$ids <- NULL else rec$ids <- target[j]
        }
        new_patches[[target[i]]] <- rec
      }
      if(!is.null(from@global[[vars[i]]])) new_global[[target[i]]] <- from@global[[vars[i]]]
    }
    provenance <- c(provenance, from@provenance)

  } else {
    if(length(vars) > 0){
      stop("layer names were given without 'from' to take them from.", call. = FALSE)
    }
    target <- character()
  }

  if(length(values) == 0 && length(target) == 0){
    stop("nothing to add.", call. = FALSE)
  }

  clash <- intersect(c(target, names(values)), names(obj@layers))
  if(length(clash) > 0){
    stop("layer(s) already present in 'obj': ", paste(clash, collapse = ", "),
         ". Use 'rename' to store them under a different name.", call. = FALSE)
  }

  for(nm in names(values)){
    if(length(values[[nm]]) != msk_ncells(obj)){
      stop("'", nm, "' has ", length(values[[nm]]), " values but the grid has ",
           msk_ncells(obj), " cells.", call. = FALSE)
    }
    new_layers[[nm]] <- values[[nm]]
  }

  .update_mosaik(obj,
                 layers = new_layers,
                 categories = new_categories,
                 patches = new_patches,
                 global = new_global,
                 provenance = provenance,
                 step = step)
}

# --- msk_pull ------------------------------------------------------------------

#' @rdname utils
#' @param layer [`character(1)`][character]\cr the layer to pull values from.
#'   Defaults to the first layer.
#' @return \code{msk_pull}: A numeric/integer vector of cell values.
#' @examples
#' # cell values as a plain vector, ready for base R
#' head(msk_pull(landscape, "canopy"))
#' mean(msk_pull(landscape, "canopy"))
#'
#' # without a layer name the first layer is pulled
#' head(msk_pull(landscape))
#' @importFrom checkmate assertClass assertCharacter
#' @export

msk_pull <- function(obj, layer = NULL){

  assertClass(x = obj, classes = "mosaik")

  if (is.null(layer)) layer <- names(obj@layers)[1]
  assertCharacter(x = layer, len = 1)

  v <- obj@layers[[layer]]
  if (is.null(v)) stop("Layer '", layer, "' not found.")
  if (inherits(v, "rle")) inverse.rle(v) else v
}

# --- .cell ----------------------------------------------------------------------

#' Find the cell that points lie in
#'
#' @param obj the mosaik.
#' @param x,y map coordinates of the points to look up.
#' @return the number of the cell each point lies in, as it indexes a pulled
#'   layer; \code{NA} for points outside the extent.
#' @noRd
.cell <- function(obj, x, y){

  assertClass(x = obj, classes = "mosaik")
  assertNumeric(x = x, any.missing = FALSE)
  assertNumeric(x = y, any.missing = FALSE, len = length(x))

  ext <- obj@extent
  res <- msk_res(obj)
  # cells are numbered row by row from the top-left corner; a point on the
  # right or lower edge of the extent belongs to the last column or row
  col <- pmin(floor((x - ext[1]) / res[1]) + 1, obj@dims[1])
  row <- pmin(floor((ext[4] - y) / res[2]) + 1, obj@dims[2])
  out <- (row - 1) * obj@dims[1] + col
  out[x < ext[1] | x > ext[2] | y < ext[3] | y > ext[4]] <- NA
  out
}

# --- msk_terra ----------------------------------------------------------------

#' Convert a mosaik to a SpatRaster
#'
#' Requires the \pkg{terra} package (listed in Suggests).
#' @param obj [`mosaik`][mosaik]\cr the mosaik to convert.
#' @return A \code{SpatRaster} with one layer per mosaik layer.
#' @examples
#' \dontrun{
#' m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
#'             vals = list(cover = sample(1:5, 100, replace = TRUE)))
#' r <- msk_terra(m)
#' }
#' @importFrom checkmate assertClass
#' @export

msk_terra <- function(obj){

  assertClass(x = obj, classes = "mosaik")

  if(!requireNamespace("terra", quietly = TRUE)){
    stop("Package 'terra' must be installed to use as_terra().", call. = FALSE)
  }

  if(length(obj@provenance) > 0){
    warning("converting mosaik to terra: provenance history will be lost.",
            call. = FALSE)
  }

  ext <- obj@extent
  dims <- obj@dims
  layerNames <- names(obj@layers)

  out <- terra::rast(ncols = dims[1], nrows = dims[2],
                     xmin = ext[1], xmax = ext[2],
                     ymin = ext[3], ymax = ext[4],
                     nlyrs = length(layerNames))

  for(i in seq_along(layerNames)){
    v <- msk_pull(obj, layerNames[i])
    terra::values(out[[i]]) <- v
  }
  names(out) <- layerNames

  # transfer categories if present
  for(i in seq_along(layerNames)){
    cats <- obj@categories[[layerNames[i]]]
    if(!is.null(cats) && !is.null(cats$gid) && !is.null(cats$val)){
      catDF <- data.frame(id = cats$gid, label = cats$val,
                          stringsAsFactors = FALSE)
      terra::set.cats(out, layer = i, value = catDF)
    }
  }

  # set CRS
  if(!is.na(obj@crs) && nchar(obj@crs) > 0){
    terra::crs(out) <- obj@crs
  }

  return(out)
}

# --- msk_spaghettify ----------------------------------------------------------

#' Turn vector geometries into a table of coordinates
#'
#' Take sf or terra vector geometries apart into the table of vertices that
#' \code{\link{msk_rasterise}} reads: coordinates without topology, the
#' "spaghetti" form of vector data. Requires the \pkg{sf} package (listed in
#' Suggests).
#' @param x an \code{sf} or \code{sfc} object, or a terra \code{SpatVector},
#'   holding one kind of geometry: points, lines or polygons, each possibly as
#'   multi-geometries.
#' @return A \code{data.frame} with one row per vertex and the columns
#'   \code{x}, \code{y}, \code{id} and \code{part}.
#' @details \code{id} is the row of the geometry in \code{x}, so attributes of
#'   \code{x} can be matched to it, and \code{part} numbers the pieces of one
#'   geometry: the parts of a multi-geometry and the rings of a polygon, holes
#'   included. Coordinates are taken as they are; reproject \code{x} to the
#'   coordinate system of the mosaik first.
#' @examples
#' \dontrun{
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"))
#' head(msk_spaghettify(nc))
#' }
#' @export

msk_spaghettify <- function(x){

  if(!requireNamespace("sf", quietly = TRUE)){
    stop("package 'sf' must be installed to use msk_spaghettify().",
         call. = FALSE)
  }
  if(inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
  if(inherits(x, "sf")) x <- sf::st_geometry(x)
  if(!inherits(x, "sfc")){
    stop("'x' must be an sf, sfc or SpatVector object.", call. = FALSE)
  }

  kinds <- unique(sub("^MULTI", "", as.character(sf::st_geometry_type(x))))
  if(length(kinds) != 1 || !kinds %in% c("POINT", "LINESTRING", "POLYGON")){
    stop("'x' must hold one kind of geometry (points, lines or polygons), ",
         "but holds ", paste(kinds, collapse = ", "), ".", call. = FALSE)
  }

  # st_coordinates numbers the nesting levels L1, L2, ..., the last of which
  # is the feature; the levels below it identify the piece within the feature
  if(length(unique(as.character(sf::st_geometry_type(x)))) > 1){
    x <- sf::st_cast(x, paste0("MULTI", kinds))
  }
  co <- sf::st_coordinates(x)
  lvl <- co[, grepl("^L[0-9]$", colnames(co)), drop = FALSE]
  if(ncol(lvl) == 0){
    id <- seq_len(nrow(co))
    part <- rep(1L, nrow(co))
  } else {
    id <- lvl[, ncol(lvl)]
    # a point needs no parts; a single line or polygon ring has one
    within <- if(ncol(lvl) > 1 && kinds != "POINT") {
      do.call(paste, as.data.frame(lvl[, -ncol(lvl), drop = FALSE]))
    } else rep("1", nrow(co))
    part <- stats::ave(match(within, unique(within)), id,
                       FUN = function(v) match(v, unique(v)))
  }

  data.frame(x = unname(co[, "X"]), y = unname(co[, "Y"]),
             id = as.numeric(id), part = as.numeric(part))
}
