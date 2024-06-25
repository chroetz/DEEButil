#' @export
writeJson <- function(obj, filePath) {
  jsonlite::write_json(
    obj,
    filePath,
    auto_unbox = TRUE,
    digits = NA,
    pretty = TRUE,
    null = "null",
    na = "string")
}


#' @export
readJson <- function(filePath) {
  jsonlite::read_json(filePath)
}
