#' @export
saveSolve <- function(A, b) {
  if (length(A) == 0) {
    res <- double(0)
    dim(res) <- c(0, NCOL(b))
    return(res)
  }
  if (all(A == 0)) {
    b[] <- NA_real_
    return(b)
  }
  res <- tryCatch(
    solve.default(A, b),
    error = function(cond) cond)
  if (!inherits(res, c("error", "condition"))) return(res)
  if (
    grepl(
      "system is exactly singular",
      res$message,
      fixed = TRUE) ||
    grepl(
      "system is computationally singular",
      res$message,
      fixed = TRUE)
    ) {
    # L2 Regularization (alternative would be MASS::ginv())
    diag(A) <- diag(A) + max(abs(A))*sqrt(.Machine$double.eps)
    return(solve.default(A, b))
  }

  signalCondition(res)
}
