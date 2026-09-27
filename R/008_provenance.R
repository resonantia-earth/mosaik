#' Record a provenance entry
#'
#' Build one provenance entry describing how a layer came to be, in the term
#' vocabulary of the W3C PROV ontology. Every mosaik-producing function calls
#' this and appends the result to \code{@provenance}, so the history travels
#' inside the object rather than in a log beside it.
#'
#' @details
#' The entry is a named list whose name is the function that produced the step,
#' and whose elements are PROV-O properties:
#'
#' \describe{
#'   \item{\code{wasGeneratedBy}}{the activity: a named list of verbs describing
#'     what was done. Always contains \code{withArguments} (the call's
#'     arguments); callers may add further verbs (\code{usedKernel},
#'     \code{appliedEquation}, ...) through \code{activity}.}
#'   \item{\code{wasDerivedFrom}}{the entity this one came from — the input
#'     layer name(s), or \code{NA} for a step that starts from nothing.}
#'   \item{\code{wasAssociatedWith}}{the agent: the software that ran the step,
#'     as \code{"<package> <version>"}, taken from the calling function's own
#'     namespace so that a downstream package is recorded as itself.}
#'   \item{\code{atTime}}{ISO-8601 timestamp with offset.}
#'   \item{\code{hash}}{digest of the activity, for change detection.}
#' }
#'
#' This is the same shape the \pkg{bitfield} package writes, so that one
#' serializer can later walk either. Neither package imports the other; the
#' shared thing is the vocabulary.
#'
#' @param fn [character(1)][character]\cr name of the function producing the
#'   step.
#' @param args [list][list]\cr the call's arguments, recorded under
#'   \code{withArguments}.
#' @param derivedFrom [character(.)][character]\cr name(s) of the layer(s) this
#'   step read, or \code{NULL} when it derives from nothing.
#' @param activity [list][list]\cr further named PROV verbs to merge into
#'   \code{wasGeneratedBy} alongside \code{withArguments}.
#' @param step [logical(1)][logical]\cr mark the entry as a replayable recipe
#'   step (see \code{\link{mdf}}). Default \code{FALSE}.
#' @return A length-one named list holding the entry.
#' @examples
#' msk_prov("mdf_scale", list(range = c(0, 800)), derivedFrom = "elevation")
#' @export

msk_prov <- function(fn, args, derivedFrom = NULL, activity = list(),
                     step = FALSE) {

  # the agent is whichever package's namespace the caller lives in, so that a
  # step run by a downstream package is recorded as that package and not as
  # mosaik. topenv() resolves correctly both when installed and under load_all.
  env <- topenv(parent.frame())
  pkg <- if (isNamespace(env)) getNamespaceName(env) else NULL
  agent <- if (is.null(pkg) || identical(pkg, "R_GlobalEnv")) {
    NA_character_
  } else {
    paste0(pkg, " ", utils::packageVersion(pkg))
  }

  generated <- c(list(withArguments = args), activity)

  entry <- list(wasGeneratedBy    = generated,
                wasDerivedFrom    = if (is.null(derivedFrom)) NA_character_ else derivedFrom,
                wasAssociatedWith = agent,
                atTime            = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                hash              = digest::digest(generated))

  if (step) entry$step <- TRUE

  structure(list(entry), names = fn)
}


#' Read the arguments out of a provenance entry
#'
#' The one place that knows where arguments sit inside an entry, so that the
#' consumers (\code{show}, \code{\link{msk_vis}}, \code{\link{mdf}}) do not each
#' hard-code the path.
#'
#' @param entry the inner list of one provenance entry.
#' @return the recorded argument list.
#' @noRd

.prov_args <- function(entry) entry$wasGeneratedBy$withArguments
