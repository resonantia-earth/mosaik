#' Categorise a mosaik
#'
#' Transform the values of a mosaik to a (smaller) set of categories, by value
#' (\code{breaks}, \code{n}) or by the share of cells each category should
#' cover (\code{shares}).
#' @param obj [`mosaik`]\cr the mosaik to modify.
#' @param breaks [`numeric(.)`][numeric]\cr the break points where categories
#'   should be delimited.
#' @param n [`integerish(1)`][integer]\cr number of categories of equal width.
#' @param shares [`numeric(.)`][numeric]\cr the share of cells in each category,
#'   from the lowest values to the highest; must sum to 1. Category 1 takes
#'   the lowest \code{shares[1]} of the cells, category 2 the next
#'   \code{shares[2]}, and so on.
#' @param layer [`character(1)`][character]\cr the layer in \code{obj} to use.
#'   Defaults to the first layer.
#' @param add [`character(1)`][character]\cr if \code{NULL} (default), overwrite
#'   \code{layer}; if a string, write to a new layer with that name.
#' @return A mosaik of the same dimension as \code{obj}, in which cells have
#'   the category value into which their values fall.
#' @details Exactly one of \code{breaks}, \code{n} and \code{shares} is given.
#'   \code{n} divides the value range of \code{layer} into \code{n} categories
#'   of equal width.
#'
#'   \code{shares} ranks the cells by their value and cuts the ranking at the
#'   cumulative shares, so the categories cover the given shares of the cells
#'   (not counting \code{NA}), up to rounding to whole cells. This turns any
#'   field into a map with known class proportions; to copy the proportions of
#'   a real map, pass \code{prop.table(table(values))}. Cells of equal value
#'   cannot be split between categories, so a cut that would fall among them is
#'   an error. This happens on layers with few distinct values, which are
#'   better categorised by \code{breaks}.
#' @examples
#' # three categories of equal width, and categories at chosen break points
#' m <- landscape |>
#'   mdf_categorise(n = 3, layer = "canopy", add = "equal_width") |>
#'   mdf_categorise(breaks = c(2, 10, 20), layer = "canopy",
#'                  add = "breaks")
#' msk_vis(m, .layer("canopy"), .layer("equal_width"), .layer("breaks"))
#'
#' # a texture turned into three classes covering 47, 30 and 23 % of the map
#' f <- mosaik(extent = c(0, 100, 0, 100), res = 1) |>
#'   drw_texture(type = "fbm", seed = 1, name = "field") |>
#'   mdf_categorise(shares = c(0.47, 0.30, 0.23), layer = "field",
#'                  add = "classes")
#' prop.table(table(msk_pull(f, "classes")))
#' msk_vis(f, .layer("field"), .layer("classes"))
#' @family operators to modify cell values
#' @importFrom checkmate assertClass assertNumeric assertIntegerish
#'   assertCharacter
#' @export

mdf_categorise <- function(obj = NULL,
                           breaks = NULL,
                           n = NULL,
                           shares = NULL,
                           layer = NULL,
                           add = NULL){

  step <- .step()
  if (.is_recipe(obj)) return(.update_mosaik(obj, step = step))

  # check arguments ----
  assertClass(x = obj, classes = "mosaik")
  assertNumeric(x = breaks, null.ok = TRUE)
  assertIntegerish(x = n, null.ok = TRUE)
  assertNumeric(x = shares, lower = 0, upper = 1, any.missing = FALSE,
                min.len = 1, null.ok = TRUE)
  assertCharacter(x = layer, null.ok = TRUE)
  assertCharacter(x = add, len = 1, null.ok = TRUE)
  if(sum(!is.null(breaks), !is.null(n), !is.null(shares)) != 1){
    stop("give exactly one of 'breaks', 'n' and 'shares'.", call. = FALSE)
  }
  # a sum other than 1 is most likely a typo; rescaling would hide it
  if(!is.null(shares) && abs(sum(shares) - 1) > 1e-8){
    stop("'shares' must sum to 1, but sum to ", signif(sum(shares), 4), ".",
         call. = FALSE)
  }

  # pull data ----
  if(is.null(layer)) layer <- names(obj@layers)[1]
  vals <- msk_pull(obj, layer)

  # body ----
  if(!is.null(shares)){

    # rank the cells and cut the ranking at the cumulative shares; ties take
    # the lowest rank of their group, so equal values share a category
    valid <- !is.na(vals)
    N <- sum(valid)
    cuts <- round(cumsum(shares) * N)
    rnk <- rank(vals[valid], ties.method = "min")
    temp <- rep(NA_integer_, length(vals))
    temp[valid] <- findInterval(rnk, cuts[-length(cuts)], left.open = TRUE) + 1L

    # cells of one value cannot be split between categories, and any rule
    # that split them would invent an order the data do not have
    sv <- sort(vals[valid])
    inner <- cuts[-length(cuts)]
    inner <- inner[inner > 0 & inner < N]
    split_at <- inner[sv[inner] == sv[inner + 1]]
    if(length(split_at) > 0){
      stop("the shares cannot be met, because the cut would split the cells ",
           "of value ", paste(unique(signif(sv[split_at], 6)), collapse = ", "),
           " between two categories.", call. = FALSE)
    }

    uVals <- sort(unique(temp[!is.na(temp)]))
    groupsLabel <- vapply(uVals, function(k) {
      r <- range(vals[!is.na(temp) & temp == k])
      paste0(round(r[1], 4), "-", round(r[2], 4))
    }, character(1))
    return(.update_mosaik(obj, values = temp, gid = uVals, val = groupsLabel,
                          step = step))
  }

  theRange <- range(vals, na.rm = TRUE)
  if(is.null(breaks)){
    assertIntegerish(x = n)
    breaks <- seq(theRange[1], theRange[2], length.out = n+1)
  }

  # the breaks must cover the values; the range is added only on a side the
  # breaks leave open, so a break beyond the values stays the outer bound
  breaks <- sort(unique(breaks))
  if(breaks[1] > theRange[1]) breaks <- c(theRange[1], breaks)
  if(breaks[length(breaks)] < theRange[2]) breaks <- c(breaks, theRange[2])
  temp <- findInterval(vals, breaks, rightmost.closed = TRUE)

  # build output ----
  uVals <- sort(unique(temp[!is.na(temp)]))
  groupsLabel <- paste0(round(breaks[uVals], 4), "-",
                        round(breaks[uVals + 1], 4))

  .update_mosaik(obj, values = temp, gid = uVals, val = groupsLabel, step = step)
}
