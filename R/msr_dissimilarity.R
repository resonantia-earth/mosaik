#' Contrast-weighted edge dissimilarity
#'
#' Weight the adjacency matrix by a user-supplied contrast matrix and
#' attach the result to the attribute table.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param contrast [`matrix`][matrix]\cr a symmetric matrix of contrast weights
#'   (values between 0 and 1). Rows and columns must be named with class
#'   values (as character). The diagonal should be 0 (no contrast within
#'   a class).
#' @param scale [`character(1)`][character]\cr the level to report at;
#'   \code{"class"} (contrast-weighted edge length per class, see
#'   \code{\link{msk_categories}}) or \code{"landscape"} (total
#'   contrast-weighted edge length, see \code{\link{msk_global}}).
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with dissimilarity values attached.
#' @details The adjacency matrix (double-counted) is computed, then
#'   multiplied element-wise by the contrast matrix. At class level, the
#'   per-class dissimilarity is the row sum of the contrast-weighted
#'   adjacency matrix. At landscape level, it is the total sum (divided
#'   by 2 to avoid double-counting).
#'
#'   In \code{\link{msr}()}, use \code{dissimilarity.class} or
#'   \code{dissimilarity.landscape}.
#' @examples
#' # build a contrast matrix for all classes in landscape
#' # (uniform contrast of 1 between all different classes)
#' cls <- as.character(sort(unique(msk_pull(landscape, "cover"))))
#' cmat <- matrix(1, nrow = length(cls), ncol = length(cls),
#'                dimnames = list(cls, cls))
#' diag(cmat) <- 0
#'
#' # class-level dissimilarity
#' m <- msr_dissimilarity(landscape, contrast = cmat, scale = "class")
#' msk_categories(m)$dissimilarity
#'
#' # landscape-level dissimilarity
#' m <- msr_dissimilarity(landscape, contrast = cmat, scale = "landscape")
#' msk_global(m)$dissimilarity
#'
#' # non-uniform contrast: some class pairs more dissimilar than others.
#' # the contrast matrix must stay symmetric, so each pair is set on both
#' # sides of the diagonal in one assignment.
#' cmat2 <- cmat
#' cmat2["44", "47"] <- cmat2["47", "44"] <- 0.2   # orchard and forest: similar
#' cmat2["35", "47"] <- cmat2["47", "35"] <- 0.9   # road and forest: very different
#' m <- msr_dissimilarity(landscape, contrast = cmat2, scale = "class")
#' msk_categories(m)$dissimilarity
#'
#' # msr: edge contrast index (dissimilarity / perimeter)
#' m <- msr_perimeter(m, scale = "class")
#' m <- msr(m, equation = "dissimilarity.class / perimeter.class",
#'          label = "edge_contrast")
#' msk_categories(m)$edge_contrast
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter assertMatrix
#' @export

msr_dissimilarity <- function(obj = NULL, contrast, scale = "class", layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertMatrix(x = contrast, mode = "numeric", min.rows = 2, min.cols = 2)
  assertChoice(x = scale, choices = c("class", "landscape"))
  assertCharacter(x = layer, null.ok = TRUE)

  if(is.null(rownames(contrast)) || is.null(colnames(contrast))){
    stop("'contrast' matrix must have row and column names matching class values.")
  }
  if(!isSymmetric(contrast)){
    stop("'contrast' matrix must be symmetric.")
  }

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)
  dims <- obj@dims

  uVals <- sort(unique(vals[!is.na(vals)]))
  uChars <- as.character(uVals)

  # check that contrast covers the classes present
  missing <- setdiff(uChars, rownames(contrast))
  if(length(missing) > 0){
    stop(sprintf("Contrast matrix is missing classes: %s",
                 paste(missing, collapse = ", ")))
  }

  # compute adjacency matrix (double-counted)
  adj <- countCellAdjacenciesCpp(vals = vals, nrow = dims[1], ncol = dims[2],
                                  doublecount = TRUE)
  rownames(adj) <- uChars
  colnames(adj) <- uChars

  # subset contrast matrix to match classes present
  cmat <- contrast[uChars, uChars]

  # element-wise product
  weighted <- adj * cmat

  if(scale == "class"){

    # per-class: row sum of weighted adjacency
    dissim <- rowSums(weighted)
    newGids <- as.integer(uVals)
    existing <- obj@categories[[layer]]
    if(!is.null(existing)){
      idx <- match(existing$gid, newGids)
      existing$dissimilarity <- ifelse(is.na(idx), NA, dissim[idx])
      obj@categories[[layer]] <- existing
    } else {
      obj@categories[[layer]] <- list(gid = newGids, dissimilarity = as.numeric(dissim))
    }

  } else {

    # landscape: total (divided by 2 because double-counted)
    obj@global[[layer]]$dissimilarity <- sum(weighted) / 2

  }

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
