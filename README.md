# SONAR: A Probabilistic Framework for Cell-Type Deconvolution in Spatial Transcriptomics

SONAR is an algorithm developed for cell-type deconvolution in spatial transcriptomics. It integrates spatial information in a balanced way to enhance performance and robustness.

## Key Features

- **Enhanced Signal-to-Noise Ratio**  
  Utilizes the similarity of spatial locations to boost the signal-to-noise ratio, helping to extract meaningful patterns from the data.

- **Local Spatial Heterogeneity Consideration**  
  Incorporates local spatial heterogeneity to avoid over-reliance on spatial information, ensuring that the inherent biological diversity is accurately captured.

- **Robust Probabilistic Framework**  
  Built on a probabilistic framework, SONAR delivers more robust and reliable results compared to deep learning–based approaches.
  
## Installation

```r
# install.packages("devtools")
devtools::install_github("lzygenomics/SONAR_R")
```

The installed package name is still `SONAR`:

```r
library(SONAR)
```

## Dependencies

SONAR_R supports two installation modes.

### Option 1: R-only native backend

This is the recommended setup for users who do not have MATLAB. It uses the new R/Rcpp implementation and runs deconvolution with `backend = "native"`. This mode does not require MATLAB, `matlabr`, or `R.matlab`.

Required:

- R version >= 4.0.5
- R packages imported by SONAR_R:
  - `this.path`
  - `Matrix`
  - `data.table`
  - `Seurat`
  - `Rcpp`
  - `gtools`
  - `ggplot2`
  - `scatterpie`
  - `viridis`
  - `corrplot`
  - `ggcorrplot`
  - `quadprog`
  - `readr`
  - `tibble`
  - `methods`
  - `reshape2`

For rendering the example R Markdown file, also install:

```r
install.packages(c("rmarkdown", "knitr"))
```

In R-only mode, the original MATLAB backend is not available. Results are written to `SONAR_results.csv` when `R.matlab` is not installed.

Example:

```r
SONAR.deconvolute(
  fname = file.path(code_path, "SONAR_main.m"),
  path = code_path,
  h = h,
  backend = "native",
  cores = 8
)
```

### Option 2: R + MATLAB compatible setup

Use this setup only if you want to keep both the new native R/Rcpp backend and the original MATLAB backend.

Additional requirements:

- MATLAB >= R2019a
- MATLAB Optimization Toolbox, required by the original `fmincon` workflow
- R packages:
  - `matlabr`
  - `R.matlab`

Install optional R packages:

```r
install.packages(c("matlabr", "R.matlab", "rmarkdown", "knitr"))
```

Native backend:

```r
SONAR.deconvolute(..., backend = "native")
```

Original MATLAB backend:

```r
SONAR.deconvolute(..., backend = "matlab")
```

## Run SONAR

1. Install the dependencies using one of the two configurations above.

2. Install SONAR_R from this repository.

3. Download or clone the SONAR_R files. This file structure helps you run the example and prepare custom data.

4. Open **Example/SONAR-entrance.Rmd** and run the example. The example uses `backend = "native"` by default and does not require MATLAB.

For running the custom dataset, you could modify the data preparation stage in **Example/SONAR-entrance.Rmd** with the same format, 
and substitute the single-cell data and spatial data in this file structure.

## Files Annotation

- **/inst/extdata/**	

  *These files store the required input information, See SONAR-entrance.Rmd for specific format requirements.*


- **/result/**

  *These files store the output results*.

  Proportions of all cell types in each spot (**SONAR.results.txt**);

  Spatial pie plots (**pie.pdf**);

  Spatial distribution of specific cell types proportion (**abs_prop.pdf** / **scaled_prop.pdf**);

  Colocalization(correlation along the spatial) for pairs of cell types (**colocalization.pdf**).


- **/core-code/**	

  *no need to operate.*

  *These files store the core code, and store the preprocessed data that delivered to SONAR*.
  
## A Brief Example

Please Follow the Example/SONAR-entrance.Rmd

## Publication Link

  https://www.nature.com/articles/s41467-023-40458-9

## How to cite SONAR

Liu, Z., Wu, D., Zhai, W. et al. SONAR enables cell type deconvolution with spatially weighted Poisson-Gamma model for spatial transcriptomics. Nat Commun 14, 4727 (2023). https://doi.org/10.1038/s41467-023-40458-9

## Thank you, and happy researching!
