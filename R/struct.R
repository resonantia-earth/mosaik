#' The \code{struct} class (S4) and its methods
#'
#' A structuring element is any 2D array with an odd number of cells in each
#' dimension and a focal cell in the middle. Each cell is compared against a
#' mosaik to probe the local pattern. Typically a conclusion is drawn based on
#' how the structuring element matches the local pattern and this conclusion is
#' stored in the value of the focal cell.
#'
#' The cells of a structuring element can have four values: NA (indifferent), 0,
#' 1, and any value > 1.
#' @slot pattern [`matrix(.)`][matrix]\cr the pattern matrix.
#' @slot rotate [`logical(1)`][logical]\cr whether the pattern shall be applied
#'   in all rotations.
#' @slot background [`numeric(1)`][numeric]\cr the value of background cells.

struct <- setClass(Class = "struct",
                   slots = c(pattern = "matrix",
                             rotate = "logical",
                             background = "numeric"
                   )
)

#' Print the structuring element
#'
#' @param object [`struct`][struct]\cr the structuring element to print.
#' @importFrom crayon yellow

setMethod(f = "show",
          signature = "struct",
          definition = function(object){

            ptrn <- object@pattern
            bg <- object@background
            rot <- object@rotate
            dims <- dim(ptrn)

            df <- data.frame(id = c(0, 1, bg), val = c("\u25ef", "\u25c9", "\u25cc"), stringsAsFactors = FALSE)

            cat(paste0("background: ", bg, "; rotate: ", ifelse(rot, "yes", "no"), "\n\n"))

            for(i in 1:dims[1]){
              temp <- as.vector(ptrn[i,])
              temp <- df$val[match(x = temp, df$id)]
              cat("   ", temp, "\n")
            }

          }
)

#' Make a structuring element
#'
#' @param type [`character(1)`][character]\cr the shape, one of 'disc', 'box'
#'   ('rectangle' and 'square' are accepted as synonyms), 'diamond' or 'cross'.
#'   Ignored when \code{custom} is given.
#' @param width [`integerish(1)`][integer]\cr width in number of cells, must be
#'   odd so that there is a focal cell.
#' @param height [`integerish(1)`][integer]\cr height in number of cells, must be
#'   odd so that there is a focal cell.
#' @param rotate [`logical(1)`][logical]\cr whether to apply all rotations.
#' @param background [`integerish(1)`][integer]\cr background cell value.
#' @param custom [`matrix(.)`][matrix]\cr a custom pattern matrix.
#' @details The generated shapes are built around the focal cell in the middle
#'   and follow the definition used by \code{mmand::shapeKernel()}. 'disc',
#'   'box' and 'diamond' are the balls of the euclidean, maximum and manhattan
#'   norm; 'cross' is not a norm ball and contains the focal row and column.
#'
#'   A cell is not treated as a point but as a square with extent, and it is
#'   included when at least half of its area falls inside the shape. This is why
#'   a 3 by 3 'disc' is filled: at that size every corner cell is more than half
#'   covered by the circle. The plus-shape is the 3 by 3 'diamond'. An unequal
#'   width and height stretch the shape along the wider axis.
#' @importFrom checkmate assertChoice assertIntegerish assertFlag assertMatrix
#' @importFrom methods new
#' @export

msk_struct <- function(type = "disc", width = 3, height = 3, rotate = FALSE,
                       background = NA, custom = NULL){

  assertFlag(x = rotate)
  assertMatrix(x = custom, null.ok = TRUE)

  if(is.na(background)) background <- NA_integer_

  if(is.null(custom)){

    assertChoice(x = type, choices = c("disc", "box", "rectangle", "square",
                                       "diamond", "cross"))
    assertIntegerish(x = width, lower = 1, len = 1, any.missing = FALSE)
    assertIntegerish(x = height, lower = 1, len = 1, any.missing = FALSE)
    if(width %% 2 == 0 | height %% 2 == 0){
      stop("'width' and 'height' must be odd, so that there is a focal cell.")
    }

    # offsets of the cell centres from the focal cell
    dx <- matrix(rep(seq_len(width) - (width + 1) / 2, each = height),
                 nrow = height, ncol = width)
    dy <- matrix(rep(seq_len(height) - (height + 1) / 2, times = width),
                 nrow = height, ncol = width)

    if(type == "cross"){

      temp <- (dx == 0 | dy == 0) * 1

    } else {

      # the shorter axis is scaled up, so the shape spans the full width and
      # height rather than only the largest ball that fits into both
      scaleX <- max(width, height) / width
      scaleY <- max(width, height) / height

      # a cell is a square, so it reaches from its near to its far edge; how
      # much of it lies inside decides whether it is included
      near <- function(d, s) s * (d - 0.5 * sign(d))
      far  <- function(d, s) s * (d + 0.5 * sign(d))

      norm <- switch(type,
                     disc = function(a, b) sqrt(a^2 + b^2),
                     diamond = function(a, b) abs(a) + abs(b),
                     function(a, b) pmax(abs(a), abs(b)))   # box

      minNorm <- norm(near(dx, scaleX), near(dy, scaleY))
      maxNorm <- norm(far(dx, scaleX), far(dy, scaleY))

      radius <- max(width, height) / 2
      coverage <- matrix(0, nrow = height, ncol = width)
      inside <- minNorm < radius
      coverage[inside] <- pmin(1, (radius - minNorm[inside]) /
                                 (maxNorm[inside] - minNorm[inside]))

      temp <- (coverage >= 0.5) * 1

    }
    dim(temp) <- c(height, width)

  } else {
    temp <- custom
  }

  out <- new(Class = "struct",
             pattern = temp,
             rotate = rotate,
             background = background)

  return(out)

}
