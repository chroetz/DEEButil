#' @export
hasValue <- function(x) {
  if (length(x) == 0) return(FALSE)
  if (length(x) > 1) return(TRUE)
  if (is.na(x)) return(FALSE)
  if (is.character(x) && nchar(x) == 0) return(FALSE)
  return(TRUE)
}
