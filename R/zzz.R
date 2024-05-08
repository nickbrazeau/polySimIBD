#------------------------------------------------
#' @title Structured Wright Fisher Simulator for Malarai Genetics
#'
#' @description Unique features of the malaria life-cycle and transmission
#'   dynamics requires extensions of typical population genetic simulators.
#'   Using a discrete-loci, discrete-time structured Wright Fisher framework, we
#'   simulate malaria population genetics forwards in time. Users are then able
#'   to capture the full Ancestral Recombination Graph.
#'
#' @docType package
#' @name polySimIBD
NULL

#------------------------------------------------
#' @useDynLib polySimIBD, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("polySimIBD", libpath)  # nocov
}
