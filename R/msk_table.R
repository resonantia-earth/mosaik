#' The values measured on a layer
#'
#' Return the class table of a layer together with its values for the whole
#' layer: the class codes, labels and colours, every value measured per class
#' by the \code{msr_*} functions and \code{\link{msr}}, and every value
#' \code{msr()} stored for the whole layer.
#' @param obj [`mosaik`]\cr the mosaik to read.
#' @param layer [`character(1)`][character]\cr the layer whose values to
#'   return. Defaults to the first layer.
#' @return A list of class \code{msk_table}. It has a print method that shows
#'   the values per class as a table, the class by class matrices by their
#'   size and the values for the whole layer below, but it is a list, so every
#'   value is reached by its name with \code{$} or \code{[[}: \code{t$area} is
#'   one value per class, in the order of \code{t$gid}, \code{t$distance} a
#'   class by class matrix, and \code{t$lpi} the one value of the whole layer.
#'   The names are those \code{\link{msr}} reads.
#' @examples
#' m <- landscape |>
#'   msr_area(layer = "cover") |>
#'   msr(equation = "area.self / sum(area.all) * 100", label = "pland",
#'       layer = "cover") |>
#'   msr(equation = "-sum(pland.all / 100 * log(pland.all / 100))",
#'       label = "shannon", layer = "cover")
#' t <- msk_table(m, layer = "cover")
#' t
#'
#' # it is a list
#' t$pland
#' t$shannon
#' @family utilities
#' @importFrom checkmate assertClass assertCharacter
#' @export

msk_table <- function(obj, layer = NULL){

  assertClass(x = obj, classes = "mosaik")
  assertCharacter(x = layer, len = 1, null.ok = TRUE)

  if (is.null(layer)) layer <- names(obj@layers)[1]
  if (!layer %in% names(obj@layers)) {
    stop("layer '", layer, "' not found in 'obj'.", call. = FALSE)
  }

  out <- c(obj@categories[[layer]], obj@global[[layer]])
  if (is.null(out)) out <- list()
  structure(out, class = c("msk_table", "list"), layer = layer,
            overall = names(obj@global[[layer]]))
}


#' @export

print.msk_table <- function(x, ...){

  overall <- attr(x, "overall")
  tbl <- unclass(x)
  attributes(tbl) <- list(names = names(x))
  cat("values of layer '", attr(x, "layer"), "'\n", sep = "")
  if (!length(tbl)) {
    cat("  --\n")
    return(invisible(x))
  }

  # one value per class goes into the table, in the order of the classes
  perClass <- setdiff(names(tbl), overall)
  n <- length(tbl$gid)
  column <- vapply(tbl[perClass], function(v) {
    n > 0 && is.atomic(v) && !is.matrix(v) && length(v) == n
  }, logical(1))
  if (any(column)) {
    print(as.data.frame(tbl[perClass[column]], stringsAsFactors = FALSE),
          row.names = FALSE)
  }

  # what does not fit a column is named with its size
  for (nm in perClass[!column]) {
    v <- tbl[[nm]]
    cat(nm, ": ", if (is.matrix(v)) {
      sprintf("%d x %d class by class matrix", nrow(v), ncol(v))
    } else {
      paste(class(v)[1], "of length", length(v))
    }, "\n", sep = "")
  }

  for (nm in overall) {
    v <- tbl[[nm]]
    cat("overall: ", nm, " = ", if (is.atomic(v) && length(v) == 1) {
      format(v)
    } else {
      paste(class(v)[1], "of length", length(v))
    }, "\n", sep = "")
  }

  invisible(x)
}
