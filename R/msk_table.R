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


#' @param x [`msk_table`][msk_table]\cr the table to print.
#' @param ... not used.
#' @param n [`integer(1)`][integer]\cr the number of classes printed; the rest
#'   are counted.
#' @rdname msk_table
#' @importFrom crayon yellow make_style has_color
#' @export

print.msk_table <- function(x, ..., n = 20){

  overall <- attr(x, "overall")
  tbl <- unclass(x)
  attributes(tbl) <- list(names = names(x))
  nClass <- length(tbl$gid)

  # the labels in a left column, as in the print of a mosaik
  lab <- function(s) yellow(formatC(s, width = -12))
  blank <- strrep(" ", 12)
  cat(lab("msk_table"), "layer ", attr(x, "layer"), " | ", nClass,
      if (nClass == 1) " class" else " classes", "\n", sep = "")
  if (!length(tbl)) return(invisible(x))

  # one value per class goes into the table, in the order of the classes
  perClass <- setdiff(names(tbl), c(overall, "gid", "val", "colour"))
  column <- vapply(tbl[perClass], function(v) {
    nClass > 0 && is.atomic(v) && !is.matrix(v) && length(v) == nClass
  }, logical(1))

  if (nClass > 0) {
    shown <- seq_len(min(nClass, n))
    # numbers are rounded for the print only; the list keeps them as they are
    cols <- lapply(tbl[perClass[column]], function(v) .format_column(v[shown]))
    cells <- c(list(gid = as.character(tbl$gid[shown])), cols)
    width <- vapply(names(cells), function(nm) max(nchar(c(nm, cells[[nm]]))), 1)
    right <- function(s, w) formatC(s, width = w)

    # the class: its label, after a square in its colour where the console
    # shows colour
    val <- if (is.null(tbl$val)) rep("", nClass) else as.character(tbl$val)
    val[is.na(val)] <- ""
    named <- any(nzchar(val))   # a layer of patches has no labels
    swatch <- named && has_color() && !is.null(tbl$colour)
    classW <- max(nchar(c("class", val[shown]))) + if (swatch) 2 else 0
    classCell <- function(i) {
      if (!named) return(NULL)
      text <- formatC(val[i], width = -(classW - if (swatch) 2 else 0))
      col <- if (swatch) tbl$colour[i] else NA
      if (!is.na(col)) paste0(make_style(col)("■"), " ", text)
      else if (swatch) paste0("  ", text)
      else text
    }

    hdr <- c(right("gid", width[["gid"]]),
             if (named) formatC("class", width = -classW),
             vapply(names(cols), function(nm) right(nm, width[[nm]]), ""))
    cat(lab("per class"), paste(hdr, collapse = "  "), "\n", sep = "")
    for (i in shown) {
      row <- c(right(cells$gid[i], width[["gid"]]), classCell(i),
               vapply(names(cols), function(nm) right(cols[[nm]][i], width[[nm]]), ""))
      cat(blank, paste(row, collapse = "  "), "\n", sep = "")
    }
    if (nClass > length(shown)) {
      rest <- nClass - length(shown)
      cat(blank, "… ", rest, if (rest == 1) " more class" else " more classes",
          "\n", sep = "")
    }
  }

  # what does not fit a column is named with its size
  rest <- perClass[!column]
  mats <- rest[vapply(tbl[rest], is.matrix, logical(1))]
  if (length(mats)) {
    cat(lab("matrices"), paste0(mats, " (", vapply(tbl[mats], function(v)
      paste(dim(v), collapse = " x "), ""), ")", collapse = ", "), "\n", sep = "")
  }
  other <- setdiff(rest, mats)
  if (length(other)) {
    cat(lab("other"), paste0(other, " (", vapply(tbl[other], function(v)
      paste(class(v)[1], "of length", length(v)), ""), ")", collapse = ", "),
      "\n", sep = "")
  }

  # the values for the whole layer, one per line
  for (i in seq_along(overall)) {
    v <- tbl[[overall[i]]]
    shownV <- if (is.atomic(v) && length(v) == 1) {
      .format_column(v)
    } else {
      paste(class(v)[1], "of length", length(v))
    }
    cat(if (i == 1) lab("overall") else blank, overall[i], " = ", shownV, "\n",
        sep = "")
  }

  invisible(x)
}

# A column of values as text. Numbers get the same number of decimals, enough
# to show about four digits of the largest one, so the decimal points line up:
# 1.67 and 24.73, 5600 and 83100.
.format_column <- function(v) {
  if (!is.numeric(v)) return(as.character(v))
  fin <- abs(v[is.finite(v)])
  if (!length(fin) || all(fin == round(fin))) {
    digits <- 0
  } else {
    digits <- min(max(3 - floor(log10(max(fin))), 0), 6)
  }
  out <- formatC(v, format = "f", digits = digits, big.mark = "")
  out[is.na(v)] <- "NA"
  out[is.infinite(v)] <- ifelse(v[is.infinite(v)] > 0, "Inf", "-Inf")
  out
}
