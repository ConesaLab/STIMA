#' selectCoord(image)
#'
#' Allows the user to select multiple points on the provided image by clicking on it.
#'
#' @param image Image where landmarks are selected 
#' @return list containing the coordinates of the selected points with two elements: x and y.
#' @export
selectCoord <- function(image) {
  plot(image)
  coordinates <- graphics::locator(type = "p")
  return(coordinates)
}


#'
#'
download_example_data <- function(dest_dir = tempdir()) {
  zip_url <- "https://github.com/vagm110901/STIMA-data/raw/main/inst/extdata/"
  destfile <- file.path(dest_dir, "exampleData.zip")
  
  utils::download.file(zip_url, destfile, mode = "wb")
  
  utils::unzip(destfile, exdir = dest_dir)
  
  return(dest_dir)
}
