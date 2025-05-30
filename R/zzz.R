.onLoad <- function(libname, pkgname) {
  if (utils::packageVersion("Seurat") != "5.0.2" | utils::packageVersion("SeuratObject") != "5.0.2") {
    utils::remove.packages("Seurat")
    utils::remove.packages("SeuratObject")
    remotes::install_version("Seurat", version = "5.0.2")
    remotes::install_version("SeuratObject", version = "5.0.2")
  }
}
