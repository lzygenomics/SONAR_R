# SONAR: A Probabilistic Framework for Cell-Type Deconvolution in Spatial Transcriptomics

SONAR is an algorithm developed for cell-type deconvolution in spatial transcriptomics. It integrates spatial information in a balanced way to enhance performance and robustness.

## Key Features

- **Enhanced Signal-to-Noise Ratio**  
  Utilizes the similarity of spatial locations to boost the signal-to-noise ratio, helping to extract meaningful patterns from the data.

- **Local Spatial Heterogeneity Consideration**  
  Incorporates local spatial heterogeneity to avoid over-reliance on spatial information, ensuring that the inherent biological diversity is accurately captured.

- **Robust Probabilistic Framework**  
  Built on a probabilistic framework, SONAR delivers more robust and reliable results compared to deep learning–based approaches.
  
## Setup 1: R-only Native Backend

Use this setup if you do not have MATLAB. This is the recommended route for `SONAR_R`.

The native backend does not call MATLAB, `fmincon`, or any `.m` file. Deconvolution is run inside R with R/Rcpp code. The installed package name is still `SONAR`.

### Step 1. Create a clean conda environment

```bash
conda create -y -n SONAR_R -c conda-forge -c bioconda \
  r-base=4.4 r-remotes r-rmarkdown r-knitr r-here r-seurat r-rcpp \
  r-this.path r-data.table r-gtools r-ggplot2 r-scatterpie r-viridis \
  r-corrplot r-ggcorrplot r-quadprog r-readr r-tibble r-reshape2 \
  pandoc git cxx-compiler c-compiler fortran-compiler make \
  pkg-config fontconfig freetype
```

### Step 2. Activate the environment

```bash
conda activate SONAR_R
```

### Step 3. Install SONAR_R

```bash
Rscript -e 'remotes::install_github("lzygenomics/SONAR_R", dependencies = FALSE, upgrade = "never")'
```

### Step 4. Verify the installation

```bash
Rscript -e 'library(SONAR); packageVersion("SONAR"); stopifnot(exists("SONAR.deconvolute.native"))'
```

This should load `SONAR` and confirm that the native backend function is available.

### Step 5. Run the example R Markdown file

```bash
git clone https://github.com/lzygenomics/SONAR_R.git
cd SONAR_R
Rscript -e 'rmarkdown::render("Example/SONAR-entrance.Rmd", output_file = "SONAR-entrance-native.html")'
```

Expected outputs:

- `Example/SONAR-entrance-native.html`
- `core-code/SONAR_results.csv`
- `result/SONAR.results.txt`
- `result/pie.pdf`, `result/abs_prop.pdf`, `result/scaled_prop.pdf`, `result/colocalization.pdf`

### Step 6. Use the native deconvolution backend

For new R scripts, call the native backend directly:

```r
SONAR.deconvolute.native(
  path = code_path,
  h = h,
  cores = 8
)
```

`cores` controls spot-wise parallel optimization. If not set, SONAR_R uses one fewer than the detected number of physical CPU cores on non-Windows systems.

For backward compatibility, this also works:

```r
SONAR.deconvolute(
  fname = "",
  path = code_path,
  h = h,
  backend = "native",
  cores = 8
)
```

When `backend = "native"`, the `fname` argument is ignored.

If you do not use conda, install R >= 4.0.5, a working C/C++ compiler, Pandoc, and all R packages listed in `DESCRIPTION`. The conda route above is the tested route.

## Setup 2: R + MATLAB Compatible Backend

Use this setup only if you want to keep both backends:

- Native R/Rcpp backend in `SONAR_R`
- Original MATLAB backend from the archived SONAR workflow

First complete **Setup 1**. Then install the optional R packages:

```bash
conda activate SONAR_R
Rscript -e 'install.packages(c("matlabr", "R.matlab"), repos = "https://cloud.r-project.org")'
```

You also need:

- MATLAB >= R2019a
- MATLAB Optimization Toolbox
- A working `matlab` command on your shell `PATH`

Native backend:

```r
SONAR.deconvolute.native(path = code_path, h = h, cores = 8)
```

Original MATLAB backend:

```r
SONAR.deconvolute(
  fname = file.path(code_path, "SONAR_main.m"),
  path = code_path,
  h = h,
  backend = "matlab",
  wait = TRUE
)
```

## Run SONAR With Your Own Data

Start from `Example/SONAR-entrance.Rmd`. Replace the single-cell expression, cell-type annotation, spatial expression, and spot coordinate files with your own data while keeping the same input format.

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
