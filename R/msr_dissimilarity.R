#' Contrast-weighted edge dissimilarity
#'
#' Weight the adjacency matrix by a user-supplied contrast matrix and
#' attach the result to the attribute table.
#' @param obj [mosaik]\cr the mosaik to measure.
#' @param contrast [matrix][matrix]\cr a symmetric matrix of contrast weights
#'   (values between 0 and 1). Rows and columns must be named with class
#'   values (as character). The diagonal should be 0 (no contrast within
#'   a class).
#' @param scale [character(1)][character]\cr the level to report at;
#'   \code{"class"} (contrast-weighted edge length per class, attached to
#'   \code{@categories}) or \code{"landscape"} (total contrast-weighted
#'   edge length, attached to \code{@global}).
#' @param layer [character(1)][character]\cr the layer to use.
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
#' # build a contrast matrix for all 9 classes in landscape
#' # (uniform contrast of 1 between all different classes)
#' cls <- c("1","11","21","24","27","31","41","44","47")
#' cmat <- matrix(1, nrow = 9, ncol = 9, dimnames = list(cls, cls))
#' diag(cmat) <- 0
#'
#' # class-level dissimilarity
#' m <- msr_dissimilarity(landscape, contrast = cmat, scale = "class")
#' m@categories$cover$dissimilarity
#'
#' # landscape-level dissimilarity
#' m <- msr_dissimilarity(landscape, contrast = cmat, scale = "landscape")
#' m@global$dissimilarity
#'
#' # non-uniform contrast: some class pairs more dissimilar than others.
#' # the contrast matrix must stay symmetric, so each pair is set on both
#' # sides of the diagonal in one assignment.
#' cmat2 <- cmat
#' cmat2["44", "47"] <- cmat2["47", "44"] <- 0.2   # similar classes
#' cmat2["1", "47"]  <- cmat2["47", "1"]  <- 0.9   # very different
#' m <- msr_dissimilarity(landscape, contrast = cmat2, scale = "class")
#' m@categories$cover$dissimilarity
#'
#' # msr: edge contrast index (dissimilarity / perimeter)
#' m <- msr_perimeter(m, scale = "class")
#' m <- msr(m, equation = "dissimilarity.class / perimeter.class",
#'          label = "edge_contrast")
#' m@categories$cover$edge_contrast
#' @family measure
#' @importFrom checkmate assertClass assertChoice assertCharacter assertMatrix
#' @export

msr_dissimilarity <- function(obj, contrast, scale = "class", layer = NULL){

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
    obj@global$dissimilarity <- sum(weighted) / 2

  }

  # provenance
  prov <- msk_prov("msr_dissimilarity", list(scale = scale, layer = layer))
  obj@provenance <- c(obj@provenance, list(prov))

  return(obj)
}
