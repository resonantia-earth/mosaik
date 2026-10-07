#' Print mosaik in the console
#'
#' @param object [`mosaik`]\cr object to \code{show}.
#' @importFrom utils head
#' @importFrom crayon yellow

setMethod(f = "show",
          signature = "mosaik",
          definition = function(object){

            layerNames <- names(object@layers)
            nLayers <- length(layerNames)
            nCells <- prod(object@dims)

            if (is.na(object@crs)) {
              myCrs <- "cartesian"
            } else {
              myCrs <- object@crs
            }

            # one line per layer with what is measured on it: per class, then
            # for the whole layer
            measured <- character()
            for (nm in layerNames) {
              perClass <- .measured_names(object@categories[[nm]])
              whole <- names(object@global[[nm]])
              parts <- c(if (length(perClass)) paste0("per class: ", paste0(perClass, collapse = ", ")),
                         if (length(whole)) paste0("overall: ", paste0(whole, collapse = ", ")))
              if (length(parts)) measured[nm] <- paste0(parts, collapse = "; ")
            }

            if (length(measured)) {
              width <- max(nchar(names(measured)))
              lines <- paste0(formatC(names(measured), width = -width), "  ", measured)
              myAttributes <- paste0(c(" ", rep("            ", length(lines) - 1)),
                                     lines, "\n")
            } else {
              myAttributes <- " --\n"
            }

            theRes <- c((object@extent[2] - object@extent[1]) / object@dims[1],
                        (object@extent[4] - object@extent[3]) / object@dims[2])

            cat(yellow("mosaik"), "\n", sep = "")
            cat("            ", nLayers, " ",
                ifelse(nLayers == 1, "layer", "layers"), " | ",
                nCells, " cells (",
                object@dims[1], "x", object@dims[2], ")\n", sep = "")
            cat(yellow("crs         "), myCrs, "\n", sep = "")
            cat(yellow("layers      "), paste0(layerNames, collapse = ", "), "\n", sep = "")
            cat(yellow("measured   "), myAttributes, sep = "")
            cat(yellow("resolution  "), theRes[1], " ", theRes[2], " (x, y)\n", sep = "")
            cat(yellow("extent     "), object@extent, "(xmin, xmax, ymin, ymax)", sep = " ")

            # show provenance if present
            steps <- .format_history(object@provenance)
            if (length(steps) > 0) {
              cat("\n")
              cat(yellow("provenance  "),
                  paste0(steps, collapse = "\n            "), "\n", sep = "")
            }
          }
)
