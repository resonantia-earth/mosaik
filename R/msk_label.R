#' Label the classes of a layer
#'
#' Give the classes of a layer a label and, optionally, a colour. Legends of
#' \code{\link{msk_vis}} show the labels and draw the classes in their colours.
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param labels [`data.frame`][data.frame]\cr one row per class, with the
#'   columns \code{id} (the class code, as in the cells of \code{layer}) and
#'   \code{label}, and optionally \code{colour}.
#' @param layer [`character(1)`][character]\cr the layer whose classes are
#'   labelled. Defaults to the first layer.
#' @return The mosaik with the labels (and colours) in the class table of
#'   \code{layer} (see \code{\link{msk_table}}).
#' @details Classes that \code{labels} does not name keep their label and
#'   colour; a class without a label so far is labelled with its code. Every
#'   \code{id} must occur in \code{layer}. Labels and colours describe the
#'   classes, so they are not measured values and survive every operator that
#'   keeps the classes.
#' @examples
#' f <- mdf_filter(landscape, cover == 47, add = "forest")
#' f <- msk_label(f, labels = data.frame(id = c(0, 1),
#'                                       label = c("open land", "forest"),
#'                                       colour = c("#e3d39a", "#1f5f2e")),
#'                layer = "forest")
#' msk_table(f, layer = "forest")
#' msk_vis(f, .layer("forest"))
#' @family utilities
#' @importFrom checkmate assertClass assertDataFrame assertCharacter
#'   assertSubset
#' @export

msk_label <- function(obj = NULL, labels, layer = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertDataFrame(x = labels, min.rows = 1)
  assertSubset(x = c("id", "label"), choices = names(labels))
  assertSubset(x = names(labels), choices = c("id", "label", "colour"))
  assertCharacter(x = layer, null.ok = TRUE)

  if (is.null(layer)) layer <- names(obj@layers)[1]
  if (!layer %in% names(obj@layers)) {
    stop("layer '", layer, "' not found in 'obj'.", call. = FALSE)
  }
  if (anyNA(labels$id) || anyDuplicated(labels$id)) {
    stop("'labels$id' must name each class once, without NA.", call. = FALSE)
  }
  vals <- msk_pull(obj, layer)
  classes <- sort(unique(vals[!is.na(vals)]))
  absent <- setdiff(labels$id, classes)
  if (length(absent)) {
    stop("class(es) ", paste(absent, collapse = ", "), " do not occur in layer '",
         layer, "'.", call. = FALSE)
  }

  # body ----
  tbl <- obj@categories[[layer]]
  if (is.null(tbl$gid)) tbl <- c(list(gid = as.integer(classes)), tbl)
  at <- match(labels$id, tbl$gid)

  if (is.null(tbl$val)) tbl$val <- as.character(tbl$gid)
  tbl$val[at] <- as.character(labels$label)
  if (!is.null(labels$colour)) {
    if (is.null(tbl$colour)) tbl$colour <- rep(NA_character_, length(tbl$gid))
    tbl$colour[at] <- as.character(labels$colour)
  }

  cats <- obj@categories
  cats[[layer]] <- tbl

  # build output ----
  .update_mosaik(obj, categories = cats, step = step)
}
