#' @export
getMonomialExponents <- function(dimension, degree) {
  degVecs <- as.matrix(expand.grid(rep(list(0:degree), dimension)))
  degVecs <- degVecs[rowSums(degVecs) <= degree, , drop=FALSE]
  return(degVecs)
}

#' @export
numberOfTermsInPoly <- \(degree, dim) {
  sum(sapply(0:degree, \(deg) choose(dim+deg-1, deg)))
}
