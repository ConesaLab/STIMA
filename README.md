
# STIMA

<!-- badges: start -->

<!-- badges: end -->

[**STIMA**](https://github.com/ConesaLab/STIMA) (*Spatial
Transcriptomics Image-based Methods for Alignment*), is designed to
align two or more ST slices or samples, enabling the comparison and
analysis of gene expression within the same regions. STIMA performs
alignment in a pairwise comparison manner, considering one slice as a
reference, which includes both the tissue microscope image and the
spatial spot matrix of gene expression. However, STIMA relies
exclusively on image data for alignment without incorporating gene
expression data, thereby preserving the independence of transcriptomic
information across samples.

STIMA includes three distinct alignment approaches:

- **Geometric Transformation Estimation Model (GTEM)**: [Krasheninnikov
  & Potapov 2012](https://doi.org/10.1134/S105466181202006X); [Zhao et
  al. 2022](https://doi.org/10.1016/B978-0-12-824349-7.00021-9)
- **Procrustes Transformation**: [Gower
  1975](https://doi.org/10.1007/BF02291478); [Murphy et
  al. 2020](https://doi.org/10.1214/19-BA1179)
- **The ImageJ plugin Register Virtual Stack Slices (RVSS-ImageJ)**:
  [ImageJ Wiki
  2022](https://imagej.github.io/plugins/register-virtual-stack-slices)

**STIMA** takes its name from the Valencian/Catalan word *estima*,
meaning *‘love’*.

## Installation

You can install STIMA from [GitHub](https://github.com/) with:

``` r
#install.packages("remotes")
remotes::install_github("ConesaLab/STIMA")
```

### ⚠️ WARNING:

If, during the installation of STIMA, R asks to update some packages to
their latest versions, proceed with caution.

STIMA has been developed to work with Seurat version 5.0.2 and not newer
versions, as significant changes were introduced in the way Spatial
Transcriptomics data is read. The semla package was used in version
1.3.1, and spacexr in version 2.2.1. If you encounter issues with any
other packages in their latest versions, consider reverting to an
earlier version to ensure STIMA functions correctly.

For development, R version 4.3.2 or earlier was used.

### Renv

If you are familiar with
[renv](https://rstudio.github.io/renv/index.html), you can install all
necessary R packages (with exact versions) with the `renv.lock` file
provided.

``` r
#install.packages("renv")
renv::restore()
```

## Documentation

### Vignettes

STIMA provides the following vignettes to better understand the
functionality of the package:

- [STIMA Usage](vignettes/STIMA_align_intrapatient.Rmd): Guide on how to
  use STIMA for aligning multiple slices. It covers three main steps: 1)
  STIMA alignment, 2) Evaluation, and 3) Deconvolution.
- [STIMA RDS-AnnData compatibility](vignettes/STIMA_PostAlignment.Rmd):
  Guide on how to process data after performing STIMA alignment. It
  covers splitting the single object that contains all alignment slices
  (the output of STIMA alignment) into separate objects—one for each
  slice—and saving CSV files for creating an AnnData object to perform
  further analyses in Python.
- [Creating Seurat objects for intra-patient
  alignment](vignettes/Rscript_Yadav2023_Merge4Slices.Rmd): Instructions
  on how to load and merge multiple 10x Visium ST slices from the same
  patient into a single Seurat object for STIMA alignment.
- [Creating Seurat objects for inter-patient
  alignment](vignettes/Rscript_Yadav2023_Merge4Patients.Rmd):
  Instructions on how to load and merge multiple 10x Visium ST samples
  from different patients into a single Seurat object for STIMA
  alignment.

### Additional useful information and tips

- It is recommended to use the **filtered H5 files** rather than the raw
  ones.
- If you crop the data and images before alignment, ensure that the
  tissue region containing gene expression information is correctly
  selected. Otherwise, the crop will be inaccurate, and the resulting
  tissue section could appear smaller than it actually is.
- After performing the alignment, you will obtain a merged RDS object
  containing both the reference image and the aligned sample (note that
  the unaligned coordinates are also retained). You can split this
  object using the Seurat function SplitObject(). Then, remove the
  reference object and save the remaining data. It is explained in
  vignette [STIMA RDS-AnnData
  compatibility](vignettes/STIMA_PostAlignment.Rmd).

### Example Data

The data used to replicate the vignettes originate from [Yadav et
al. (2023)](https://doi.org/10.1016/j.neuron.2023.01.007). To access the
data used in this package, you can visit the repository
[STIMA-exampleData](https://github.com/vagm110901/STIMA-exampleData), or
download the raw files directly from the Gene Expression Omnibus (GEO)
under the accession numbers GEO:
[GSE190442](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190442)
and GEO:
[GSE222322](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222322).
