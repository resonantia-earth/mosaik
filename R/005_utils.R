#' Utility functions that act on the tabular information of a mosaik
#'
#' These utility functions are meant to immitate the tidy dplyr logic on the
#' tables of a mosaik. \code{msk_select} keeps only the specified layers,
#' \code{msk_remove} drops them, \code{msk_add} brings layers in from a second
#' mosaik on the same grid. \code{msk_pull} extracts cell values from a layer as
#' a vector (like \code{dplyr::pull}). To select cells rather than layers, see
#' \code{\link{mdf_filter}}.
#'
#' @param obj [mosaik][mosaik]\cr the mosaik object.
#' @param ... layer names (unquoted or character strings). For
#'   \code{msk_select} and \code{msk_remove} these name layers of \code{obj},
#'   for \code{msk_add} they name layers of \code{from}.
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
#' msk_select(landscape, "cover", "intensity")
#' @importFrom checkmate assertClass
#' @export

msk_select <- function(obj, ...){

  assertClass(x = obj, classes = "mosaik")

  vars <- as.character(match.call(expand.dots = FALSE)$...)

  layerNames <- names(obj@layers)
  keep <- vars[vars %in% layerNames]
  if(length(keep) == 0){
    stop("none of the requested layers found: ", paste(vars, collapse = ", "))
  }

  new_layers <- obj@layers[keep]
  new_categories <- obj@categories[intersect(keep, names(obj@categories))]

  prov <- msk_prov("msk_select", list(layers = paste(keep, collapse = ", ")))

  new_mosaik(
    extent     = obj@extent,
    dims       = obj@dims,
    layers     = new_layers,
    categories = new_categories,
    patches    = obj@patches,
    global     = obj@global,
    crs        = obj@crs,
    provenance = c(obj@provenance, list(prov))
  )
}

# --- msk_remove ---------------------------------------------------------------

#' @rdname utils
#' @return \code{msk_remove}: A mosaik without the named layers.
#' @examples
#' # drop a layer, keeping the rest
#' msk_remove(landscape, intensity)
#'
#' # removing every layer is refused
#' try(msk_remove(landscape, cover, intensity))
#' @importFrom checkmate assertClass
#' @export

msk_remove <- function(obj, ...){

  assertClass(x = obj, classes = "mosaik")

  vars <- as.character(match.call(expand.dots = FALSE)$...)

  layerNames <- names(obj@layers)
  drop <- vars[vars %in% layerNames]
  keep <- setdiff(layerNames, drop)
  if(length(keep) == 0){
    stop("removing ", paste(drop, collapse = ", "), " would leave no layer.")
  }

  new_layers <- obj@layers[keep]
  new_categories <- obj@categories[intersect(keep, names(obj@categories))]

  prov <- msk_prov("msk_remove", list(layers = paste(drop, collapse = ", ")))

  new_mosaik(
    extent     = obj@extent,
    dims       = obj@dims,
    layers     = new_layers,
    categories = new_categories,
    patches    = obj@patches,
    global     = obj@global,
    crs        = obj@crs,
    provenance = c(obj@provenance, list(prov))
  )
}

# --- msk_add ------------------------------------------------------------------

#' @rdname utils
#' @param from [mosaik][mosaik]\cr the mosaik to take layers from. Must sit on
#'   the same grid as \code{obj}, i.e. carry the same \code{extent},
#'   \code{dims} and \code{crs}.
#' @param rename [character(.)][character]\cr optional names under which the
#'   layers should be stored in \code{obj}, in the order they are given in
#'   \code{...}. Needed when a name is already taken.
#' @return \code{msk_add}: A mosaik with the layers of \code{from} added.
#' @details \code{msk_add} is the inverse of \code{msk_select} and the way to
#'   combine two mosaiks: operators such as \code{\link{mdf_blend}} work within
#'   a single object, so a layer computed elsewhere is brought in first and
#'   combined afterwards. The grid is checked here, which is why no operator
#'   needs to check it again. Any \code{@categories} entry belonging to an added
#'   layer travels with it; \code{@patches} and \code{@global} do not, since
#'   they describe the landscape a measurement ran on rather than a single
#'   layer.
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
#' @importFrom checkmate assertClass assertCharacter
#' @export

msk_add <- function(obj, from, ..., rename = NULL){

  assertClass(x = obj, classes = "mosaik")
  assertClass(x = from, classes = "mosaik")
  assertCharacter(x = rename, null.ok = TRUE, any.missing = FALSE)

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
  vars <- as.character(match.call(expand.dots = FALSE)$...)
  if(length(vars) == 0) vars <- names(from@layers)

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

  clash <- intersect(target, names(obj@layers))
  if(length(clash) > 0){
    stop("layer(s) already present in 'obj': ", paste(clash, collapse = ", "),
         ". Use 'rename' to store them under a different name.", call. = FALSE)
  }

  # body ----
  new_layers <- obj@layers
  new_categories <- obj@categories
  for(i in seq_along(vars)){
    new_layers[[target[i]]] <- from@layers[[vars[i]]]
    cats <- from@categories[[vars[i]]]
    if(!is.null(cats)) new_categories[[target[i]]] <- cats
  }

  prov <- msk_prov("msk_add", list(layers = paste(vars, collapse = ", "),
                                     as = paste(target, collapse = ", ")))

  new_mosaik(
    extent     = obj@extent,
    dims       = obj@dims,
    layers     = new_layers,
    categories = new_categories,
    patches    = obj@patches,
    global     = obj@global,
    crs        = obj@crs,
    provenance = c(obj@provenance, from@provenance, list(prov))
  )
}

# --- msk_pull ------------------------------------------------------------------

#' @rdname utils
#' @param layer [character(1)][character]\cr the layer to pull values from.
#'   Defaults to the first layer.
#' @return \code{msk_pull}: A numeric/integer vector of cell values.
#' @examples
#' # cell values as a plain vector, ready for base R
#' head(msk_pull(landscape, "intensity"))
#' mean(msk_pull(landscape, "intensity"))
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

# --- msk_terra ----------------------------------------------------------------

#' Convert a mosaik to a SpatRaster
#'
#' Requires the \pkg{terra} package (listed in Suggests).
#' @param obj [mosaik][mosaik]\cr the mosaik to convert.
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
