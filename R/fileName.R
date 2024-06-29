#' @export
getUniqueFileName <- function(dirPath=".", prefix=NULL, fileExtension="", identifyingObject = NULL, timeStamp = FALSE, fullPath = FALSE, createDir = TRUE) {
  if (createDir && !dir.exists(dirPath)) dir.create(dirPath)
  if (timeStamp) {
    timeStampStr <- stringr::str_replace(format(Sys.time(), "%Y-%m-%d-%H-%M-%OS6"), stringr::fixed("."), "-")
  } else {
    timeStampStr <- NULL
  }
  if (is.null(identifyingObject)) {
    pattern <- paste0(c(prefix, timeStampStr), collapse="_")
    if (!hasValue(pattern)) {
      pattern <- ""
    } else {
      pattern <- paste0(pattern, "_")
    }
    fileName <-
      tempfile(
        pattern = pattern,
        tmpdir = dirPath,
        fileext = fileExtension
      ) |>
      basename()
  } else {
    fileName <- paste0(
      paste(
        c(prefix,
        timeStampStr,
        rlang::hash(identifyingObject)),
        collapse = "_"
      ),
      fileExtension)
  }
  if (fullPath) {
    result <- normalizePath(file.path(dirPath, fileName), winslash = "/", mustWork = FALSE)
  } else {
    result <- fileName
  }
  return(result)
}
