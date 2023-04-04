#' @export
evalGaussianProcess <- function(x, locations, weights, bandwidth) {
  kernelVector <- expKernelVector(x, locations, bandwidth)
  crossprod(kernelVector, weights)
}

#' @export
calculateGaussianProcessWeights <- function(locations, values, bandwidth, regulation) {
  kernelMatrix <- expKernelMatrix(locations, bandwidth, regulation)
  weights <- solve.default(kernelMatrix, values)
  return(weights)
}
