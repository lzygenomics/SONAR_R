#source("helper_preprocess.R")
#' SONAR.preprocess
#'
#' @param sc_count scRNA-seq expression matrix of reference(Details and formats are in vignettes)
#' @param sc_cell_type cell type annotation for reference(Details and formats are in vignettes)
#' @param sc_nUMI cell library size(Details and formats are in vignettes)
#' @param sp_coords spot expression coordinate(Details and formats are in vignettes)
#' @param sp_count spot expression matrix(Details and formats are in vignettes)
#' @param sp_nUMI spot library size(Details and formats are in vignettes)
#' @param cores the number of CPU threads that you want to use
#' @param type_min_cell Specifies a threshold to filter out cell types with a smaller number of cells
#' @param spot_min_UMI Specifies a threshold to filter out spots with fewer UMIs
#' @param preclus_resolution Specifies the resolution of the preclustering
#' @param exp_cutoff_reg minimum normalized gene expression for genes to be included in the regression.
#' @param fc_cutoff_reg minimum log-fold-change (across cell types) for genes to be included in the regression.
#' @param exp_cutoff_plat minimum normalized gene expression for genes to be included in the platform effect normalization step.
#' @param fc_cutoff_plat minimum log-fold-change (across cell types) for genes to be included in the platform effect normalization step.
#' @import methods Seurat utils Matrix
#' @return preprocessed data
#' @export

SONAR.preprocess<-function(sc_count,sc_cell_type,sc_nUMI,sp_coords,sp_count,sp_nUMI,cores=8,type_min_cell=0,spot_min_UMI=100,preclus_resolution=0.7,exp_cutoff_plat = 0.000125, fc_cutoff_plat = 0.5, exp_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75){
  reference <- Reference(sc_count,sc_cell_type,sc_nUMI)
  puck <- SpatialRNA(sp_coords,sp_count,sp_nUMI)
  SONAR_obj <- create.SONAR(puck, reference, max_cores = cores, CELL_MIN_INSTANCE = type_min_cell,UMI_min = spot_min_UMI,gene_cutoff = exp_cutoff_plat, fc_cutoff = fc_cutoff_plat, gene_cutoff_reg = exp_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg)
  SONAR_obj <-fitBulk(SONAR_obj)
  #u[j,t] is preprocessed type t,gene j 's expression
  u <- SONAR_obj@cell_type_info$renorm[[1]]
  #y[j,n]is spot n,gene j 's expression
  y <- as.matrix(SONAR_obj@spatialRNA@counts)
  #N[n] spot n library size
  N <- SONAR_obj@spatialRNA@nUMI
  coord<-SONAR_obj@spatialRNA@coords
  #preclustering
  print("####Part 4: ")
  print("Begin: Pre-clustering")
  yy<-CreateSeuratObject(y)
  yy<-SCTransform(yy)
  yy <- RunPCA(yy, assay = "SCT", verbose = FALSE)
  yy <- FindNeighbors(yy,reduction = "pca", dims=1:30)
  yy <- FindClusters(object = yy, verbose = T, resolution = preclus_resolution)
  yy <- RunUMAP(yy, reduction = "pca", dims = 1:10, label = T)
  #DimPlot(yy, reduction = "umap")
  precluster_label<-as.numeric(Idents(yy))
  print("End: Pre-clustering")
  p<-list(u=u,y=y,N=N,coord=coord,precluster_label=precluster_label)
  print("All preprocessing is complete!")
  return(p)
}
#' SONAR.deliver
#'
#' @param processed_data the data from preprocessed steps
#' @param path the path of delivering the preprocessed data to SONAR
#' @import utils
#' @return this part is to save the files, return the complete state
#' @export

SONAR.deliver<-function(processed_data,path){
  write.table(processed_data$u,file = paste(path,'u.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'u.txt',sep=""))) {
    stop("reference file error")
  }
  write.table(processed_data$y,file = paste(path,'y.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'y.txt',sep=""))) {
    stop("spatial file error")
  }
  write.table(processed_data$N,file = paste(path,'N.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'N.txt',sep=""))) {
    stop("spots library file error")
  }
  write.table(processed_data$coord,file=paste(path,'coord.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'coord.txt',sep=""))) {
    stop("coord file ready error")
  }
  write.table(processed_data$precluster_label,file=paste(path,'label.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'label.txt',sep=""))) {
    stop("precluster file error")
  }
  return(print("deliver complete!"))
}
#' SONAR.deconvolute.native
#'
#' @param path the path to preprocessed data
#' @param h bandwidth
#' @param cores the number of CPU threads used for spot-wise optimization
#' @param maxit maximum iterations for each spot-wise optimizer
#' @param lower lower bound used for all optimized parameters
#' @param seed optional random seed for optimizer initialization
#' @param verbose show progress
#' @param return_result return the cell-type proportion matrix instead of a status code
#' @return the state of deconvolution
#' @export

SONAR.deconvolute.native<-function(path,h,cores = max(1, parallel::detectCores(logical = FALSE) - 1, na.rm = TRUE),
                                   maxit = 1000, lower = 1e-8, seed = NULL,
                                   verbose = TRUE, return_result = FALSE)
{
  stopifnot(dir.exists(path))
  input_files <- file.path(path, c("u.txt", "y.txt", "N.txt", "coord.txt", "label.txt"))
  if (!all(file.exists(input_files))) {
    stop("missing SONAR input files: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
  }

  read_matrix <- function(file) {
    as.matrix(utils::read.csv(file, check.names = FALSE, row.names = 1))
  }

  u_raw <- read_matrix(file.path(path, "u.txt"))
  y_raw <- read_matrix(file.path(path, "y.txt"))
  N <- as.numeric(utils::read.csv(file.path(path, "N.txt"), check.names = FALSE, row.names = 1)[[1]])
  coord <- read_matrix(file.path(path, "coord.txt"))
  label <- as.numeric(utils::read.csv(file.path(path, "label.txt"), check.names = FALSE, row.names = 1)[[1]])

  if (nrow(u_raw) != nrow(y_raw)) {
    stop("u.txt and y.txt must have the same genes in rows")
  }
  if (ncol(y_raw) != length(N) || ncol(y_raw) != nrow(coord) || ncol(y_raw) != length(label)) {
    stop("spot counts, N, coord, and label dimensions do not match")
  }

  u <- rbind(Intercept = rep(1, nrow(u_raw)), t(u_raw))
  y <- t(y_raw)
  n_spots <- length(N)

  D <- as.matrix(stats::dist(coord))
  same_label <- outer(label, label, "==")
  w <- matrix(0, n_spots, n_spots)
  neighbor <- D < h & same_label
  w[neighbor] <- (1 - (D[neighbor] / h)^2)^2

  norm_y <- y
  y_sums <- rowSums(norm_y)
  if (any(y_sums <= 0)) {
    stop("all spots must have positive total expression")
  }
  norm_y <- norm_y / y_sums
  norm_lengths <- sqrt(rowSums(norm_y^2))
  cosine_sim <- tcrossprod(norm_y) / tcrossprod(norm_lengths)
  cosine_sim <- pmax(pmin(cosine_sim, 1), -1)
  cos_dist <- 1 - cosine_sim
  w <- w^(100 * cos_dist)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  xindex <- nrow(u)
  x0 <- c(stats::runif(xindex), 10)

  optimize_one <- function(i) {
    fit <- stats::optim(
      par = x0,
      fn = sonar_negloglik_cpp,
      gr = sonar_grad_cpp,
      method = "L-BFGS-B",
      lower = rep(lower, xindex + 1),
      control = list(maxit = maxit),
      u = u,
      N = N,
      w = w,
      y = y,
      spot = i
    )
    if (fit$convergence != 0 && verbose) {
      warning("optimizer did not fully converge for spot ", i, ": code ", fit$convergence,
              call. = FALSE, immediate. = TRUE)
    }
    fit$par
  }

  if (verbose) {
    message("Running SONAR native deconvolution with ", cores, " core(s)")
  }
  if (cores > 1 && .Platform$OS.type != "windows") {
    par_list <- parallel::mclapply(seq_len(n_spots), optimize_one, mc.cores = cores)
  } else {
    par_list <- lapply(seq_len(n_spots), optimize_one)
  }

  JIE <- do.call(rbind, par_list)
  JIE <- JIE[, -c(1, ncol(JIE)), drop = FALSE]
  row_sums <- rowSums(JIE)
  if (any(row_sums <= 0)) {
    stop("native optimizer returned non-positive cell-type sums")
  }
  JIE <- JIE / row_sums
  colnames(JIE) <- colnames(u_raw)
  rownames(JIE) <- rownames(coord)

  results_file <- file.path(path, "SONAR_results.mat")
  csv_results_file <- file.path(path, "SONAR_results.csv")
  tryCatch(
    {
      if (!requireNamespace("R.matlab", quietly = TRUE)) {
        stop("R.matlab is not installed")
      }
      R.matlab::writeMat(results_file, JIE = JIE)
      if (verbose) {
        message("Data successfully saved in MAT format.")
      }
    },
    error = function(e) {
      if (verbose) {
        message("Failed to save in MAT format: ", conditionMessage(e))
        message("Saving data in CSV format instead.")
      }
      utils::write.csv(JIE, csv_results_file)
    }
  )

  if (verbose) {
    message("Complete")
  }
  if (return_result) {
    return(JIE)
  }
  return(0)
}

#' SONAR.deconvolute
#'
#' @param fname the path to SONAR's MATLAB core codes
#' @param path the path to preprocessed data
#' @param h bandwidth
#' @param verbose show the command/progress
#' @param desktop Determines whether the matlab run process is displayed
#' @param splash Determines whether the matlab run process is displayed
#' @param display Determines whether the matlab run process is displayed
#' @param wait wait the SONAR's deconvolution complete
#' @param single_thread Determine whether to use only single threads
#' @param backend deconvolution backend: native R/Rcpp or MATLAB
#' @param cores the number of CPU threads used by the native backend
#' @param maxit maximum iterations for each native spot-wise optimizer
#' @param lower lower bound used by the native optimizer
#' @param seed optional random seed for native optimizer initialization
#' @return the state of deconvolution
#' @export

SONAR.deconvolute<-function (fname,path,h,verbose = TRUE, desktop = FALSE, splash = FALSE,
                             display = FALSE, wait = FALSE, single_thread = FALSE,
                             backend = c("native", "matlab"),
                             cores = max(1, parallel::detectCores(logical = FALSE) - 1, na.rm = TRUE),
                             maxit = 1000, lower = 1e-8, seed = NULL)
{
  backend <- match.arg(backend)
  if (backend == "native") {
    if (single_thread) {
      cores <- 1
    }
    return(SONAR.deconvolute.native(path = path, h = h, cores = cores, maxit = maxit,
                                    lower = lower, seed = seed, verbose = verbose))
  }

  stopifnot(file.exists(fname))
  if (!requireNamespace("matlabr", quietly = TRUE)) {
    stop("The MATLAB backend requires the matlabr package. Use backend = 'native' to avoid MATLAB.")
  }
  matcmd = matlabr::get_matlab(desktop = desktop, splash = splash,
                               display = display, wait = wait, single_thread = single_thread)
  cmd=paste0("userpath('",path,"');","thepath='",path,"';h=",h,";","SONAR_main;","exit")
  cmd = paste0(matcmd, '\"',cmd, '\"')
  if (verbose) {
    message("Command run is:")
    message(cmd)
  }
  x <- system(cmd, wait = wait)
  return(x)
}
