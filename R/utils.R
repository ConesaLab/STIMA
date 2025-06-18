#' selectCoord(image)
#'
#' Allows the user to select multiple points on the provided image by clicking on it.
#'
#' @param image Image where landmarks are selected 
#' @return list containing the coordinates of the selected points with two elements: x and y.
#' @importFrom Seurat GetImage GetTissueCoordinates
#' @export
selectCoord <- function(image) {
  plot(image)
  coordinates <- graphics::locator(type = "p")
  return(coordinates)
}

#'cropImage(object, image_name, padding = 10)
#'
#' This function crops a square region from a spatial image within a Seurat object,
#' centered around the tissue coordinates, with optional padding. The cropped image and 
#' corresponding coordinates are rescaled and updated in the object.
#'
#' @param object A Seurat object containing spatial image data in the \code{@images} slot.
#' @param image_name A character string indicating the name of the image within the Seurat object to crop.
#' @param padding An integer specifying how many pixels of padding to add around the tissue region (default is 10).
#'
#' @return A modified Seurat object with the specified image cropped and updated.
#' @export
cropImage <- function(object, image_name, padding = 10) {

  img_obj <- object@images[[image_name]]
  img_mat <- GetImage(img_obj)
  coordslow <- GetTissueCoordinates(img_obj)
  
  # Crop coordinates including a padding
  bbox_coords <- apply(coordslow, 2, range)
  x_min <- floor(bbox_coords[1, 2]) - padding
  x_max <- ceiling(bbox_coords[2, 2]) + padding
  y_min <- floor(bbox_coords[1, 1]) - padding
  y_max <- ceiling(bbox_coords[2, 1]) + padding
  
  # square
  width <- x_max - x_min; height <- y_max - y_min
  size <- max(width, height)
  
  # center the square
  x_center <- round((x_min + x_max) / 2)
  y_center <- round((y_min + y_max) / 2)
  half_size <- round(size / 2)
  
  x_min_sq <- max(1, x_center - half_size)
  x_max_sq <- x_min_sq + size
  y_min_sq <- max(1, y_center - half_size)
  y_max_sq <- y_min_sq + size
  
  # The size must be smaller or equal than the original image 
  img_dims <- dim(img_obj@image)
  x_max_sq <- min(x_max_sq, img_dims[2])
  y_max_sq <- min(y_max_sq, img_dims[1])
  x_min_sq <- x_max_sq - size
  y_min_sq <- y_max_sq - size
  
  # Image crop
  cropped_image <- img_obj@image[y_min_sq:y_max_sq, x_min_sq:x_max_sq, , drop = FALSE]
  
  # Coordinates crop
  sf <- img_obj@scale.factors
  coordslow$imagerow <- coordslow$imagerow - y_min_sq
  coordslow$imagecol <- coordslow$imagecol - x_min_sq
  coordsfull <- coordslow / sf$lowres
  
  orig_dims <- dim(img_obj)
  crop_dims <- dim(cropped_image)
  scale_factor_ratio_y <- orig_dims[1] / crop_dims[1]

  # Save the cropped object
  new_img_obj <- img_obj
  new_img_obj@image <- cropped_image
  new_img_obj@coordinates$imagerow <- coordsfull$imagerow
  new_img_obj@coordinates$imagecol <- coordsfull$imagecol
  new_img_obj@spot.radius <- img_obj@spot.radius * scale_factor_ratio_y
  
  object@images[[image_name]] <- new_img_obj
  
  return(object)
}



#' download_example_data(dest_dir = tempdir())
#'
#' Downloads a ZIP file containing example data for the STIMA package
#' from the official GitHub repository and extracts it to the specified directory.
#'
#' @param dest_dir A character string indicating the destination directory where the
#' data should be saved and unzipped. Defaults to a temporary directory (\code{tempdir()}).
#' @return A character string with the path to the directory where the data was extracted.
#' @export
download_example_data <- function(dest_dir = tempdir()) {
  zip_url <- "https://github.com/vagm110901/STIMA-data/raw/main/inst/extdata/"
  destfile <- file.path(dest_dir, "exampleData.zip")
  
  utils::download.file(zip_url, destfile, mode = "wb")
  
  utils::unzip(destfile, exdir = dest_dir)
  
  return(dest_dir)
}
