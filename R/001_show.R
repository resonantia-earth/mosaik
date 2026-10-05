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

            # show patch and landscape results per layer
            for (s in c("patches", "global")) {
              for (nm in names(methods::slot(object, s))) {
                entry <- methods::slot(object, s)[[nm]]
                fields <- names(entry)
                what <- nm
                if (s == "patches") {
                  fields <- setdiff(fields, c("class", "patch", "ids"))
                  # how many patches, and where they are numbered
                  what <- paste0(nm, ", ", length(entry$patch), " patches",
                                 if (!is.null(entry$ids)) paste0(" in '", entry$ids, "'"))
                }
                prefix <- if (is.null(myAttributes)) " " else "            "
                myAttributes <- c(myAttributes, paste0(prefix, "(", s, ") ", what,
                                                       if (length(fields)) paste0(" (", paste0(fields, collapse = ", "), ")"),
                                                       "\n"))
              }
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
            steps <- .format_history(object@provenance)
            if (length(steps) > 0) {
              cat("\n")
              cat(yellow("provenance  "),
                  paste0(steps, collapse = "\n            "), "\n", sep = "")
            }
          }
)
