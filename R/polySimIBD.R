#------------------------------------------------
#' @title An example package
#'
#' @description An example package to demonstrate development pipeline.
#'
#' @docType package
#' @name polySimIBD
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib polySimIBD, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("polySimIBD", libpath)
}
