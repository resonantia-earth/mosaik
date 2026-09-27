#' Print mosaik in the console
#'
#' @param object [mosaik]\cr object to \code{show}.
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

            myAttributes <- NULL

            # show categories info per layer
            if (length(object@categories) > 0) {
              for (nm in names(object@categories)) {
                entry <- object@categories[[nm]]
                metricCols <- setdiff(names(entry), c("gid", "val"))
                nCats <- length(entry$gid)
                catStr <- paste0(nm, " (", nCats, " classes")
                if (length(metricCols) > 0) {
                  catStr <- paste0(catStr, "; ", paste0(metricCols, collapse = ", "))
                }
                catStr <- paste0(catStr, ")")
                if (is.null(myAttributes)) {
                  myAttributes <- c(myAttributes, paste0(" (categories) ", catStr, "\n"))
                } else {
                  myAttributes <- c(myAttributes, paste0("            (categories) ", catStr, "\n"))
                }
              }
            }

            # show patches if populated
            if (length(object@patches) > 0) {
              patchNames <- names(object@patches)
              prefix <- if (is.null(myAttributes)) " " else "            "
              myAttributes <- c(myAttributes, paste0(prefix, "(patches) ",
                                                     paste0(patchNames, collapse = ", "), "\n"))
            }

            # show global if populated
            if (length(object@global) > 0) {
              globalNames <- names(object@global)
              prefix <- if (is.null(myAttributes)) " " else "            "
              myAttributes <- c(myAttributes, paste0(prefix, "(global) ",
                                                     paste0(globalNames, collapse = ", "), "\n"))
            }

            if (is.null(myAttributes)) {
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
            cat(yellow("data       "), myAttributes, sep = "")
            cat(yellow("resolution  "), theRes[1], " ", theRes[2], " (x, y)\n", sep = "")
            cat(yellow("extent     "), object@extent, "(xmin, xmax, ymin, ymax)", sep = " ")

            # show provenance if present
            theHist <- object@provenance
            if (length(theHist) > 0) {
              cat("\n")
              steps <- vapply(theHist, function(step){
                if (is.character(step)) {
                  step
                } else {
                  fn <- names(step)[1]
                  info <- step[[fn]]
                  args <- .prov_args(info)
                  args_str <- paste0(names(args), "=", args, collapse = ", ")
                  paste0(fn, "(", args_str, ") [", info$atTime, "]")
                }
              }, character(1))
              cat(yellow("provenance  "), paste0(steps, collapse = " -> "), sep = "")
            }
          }
)
