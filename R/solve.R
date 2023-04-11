#' @export
saveSolve <- function(A, b) {
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
