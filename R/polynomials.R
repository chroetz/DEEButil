#' @export
numberOfTermsInPoly <- \(degree, dim) {
  choose(degree + dim, dim)
}
