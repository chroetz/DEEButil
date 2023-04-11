#' @export
getMonomialExponents <- function(dimension, degree) {
  degVecs <- as.matrix(expand.grid(rep(list(0:degree), dimension)))
  degVecs <- degVecs[rowSums(degVecs) <= degree, , drop=FALSE]
  return(degVecs)
}

numberOfTermsInPoly <- \(polyDeg, d) {
  sum(sapply(0:polyDeg, \(deg) choose(d+deg-1, deg)))
}
