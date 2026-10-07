#' Contrast-weighted edge dissimilarity
#'
#' Weight the adjacencies of each class by a user-supplied contrast matrix and
#' store the result in the class table of the layer.
#' @param obj [`mosaik`]\cr the mosaik to measure.
#' @param contrast [`matrix`][matrix]\cr a symmetric matrix of contrast weights
#'   (values between 0 and 1). Rows and columns must be named with class
#'   values (as character). The diagonal should be 0 (no contrast within
#'   a class).
#' @param layer [`character(1)`][character]\cr the layer to use.
#'   Defaults to the first layer.
#' @return The input mosaik with \code{dissimilarity} added to the class table
#'   of \code{layer} (see \code{\link{msk_table}}).
#' @details The adjacency matrix (double-counted) is computed, then
#'   multiplied element-wise by the contrast matrix. The dissimilarity of a
#'   class is the row sum of the contrast-weighted adjacency matrix. The total
#'   contrast-weighted edge length of the layer is
#'   \code{sum(dissimilarity.all) / 2} in \code{\link{msr}}, since every edge
#'   is counted from both sides.
#' @examples
#' # build a contrast matrix for all classes in landscape
#' # (uniform contrast of 1 between all different classes)
#' cls <- as.character(sort(unique(msk_pull(landscape, "cover"))))
#' cmat <- matrix(1, nrow = length(cls), ncol = length(cls),
#'                dimnames = list(cls, cls))
#' diag(cmat) <- 0
#'
#' m <- msr_dissimilarity(landscape, contrast = cmat)
#' msk_table(m)$dissimilarity
#'
#' # non-uniform contrast: some class pairs more dissimilar than others.
#' # the contrast matrix must stay symmetric, so each pair is set on both
#' # sides of the diagonal in one assignment.
#' cmat2 <- cmat
#' cmat2["44", "47"] <- cmat2["47", "44"] <- 0.2   # orchard and forest: similar
#' cmat2["35", "47"] <- cmat2["47", "35"] <- 0.9   # road and forest: very different
#' m <- msr_dissimilarity(landscape, contrast = cmat2)
#' msk_table(m)$dissimilarity
#'
#' # msr: edge contrast index (dissimilarity / perimeter)
#' m <- msr_perimeter(m)
#' m <- msr(m, equation = "dissimilarity.self / perimeter.self",
#'          label = "contrast")
#' msk_table(m)$contrast
#' @family measure
#' @importFrom checkmate assertClass assertCharacter assertMatrix
#' @export

msr_dissimilarity <- function(obj = NULL, contrast, layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  assertClass(x = obj, classes = "mosaik")
  assertMatrix(x = contrast, mode = "numeric", min.rows = 2, min.cols = 2)
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
  adj <- countCellAdjacenciesCpp(vals = vals, nrow = dims[2], ncol = dims[1],
                                  doublecount = TRUE)
  rownames(adj) <- uChars
  colnames(adj) <- uChars

  # subset contrast matrix to match classes present
  cmat <- contrast[uChars, uChars]

  # per class: row sum of the contrast-weighted adjacencies
  dissim <- as.numeric(rowSums(adj * cmat))
  obj <- .store_class(obj, layer, "dissimilarity", as.integer(uVals), dissim)

  # provenance
  obj <- .update_mosaik(obj, step = step)

  return(obj)
}
