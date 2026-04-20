source("reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R")

# Helper function on obtaining downsampled UMI
downsample_and_count_UMI <- function(data, ratios, target_genes = NULL, seed = 42) {
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Obtain the observed reads vector
  cell_num <- length(unique(data$cell_id))
  num_observed_reads <- sum(data$num_reads)
  reads_vec <- rep(seq_along(data$num_reads), data$num_reads)

  # Perform downsampling for all ratios
  num_downsampled_reads <- round(num_observed_reads * ratios)

  if (is.null(target_genes)) {
    # Standard downsampling - count all UMIs
    num_downsampled_UMIs <- sapply(num_downsampled_reads,
                                   function(reads) length(unique(sample(reads_vec, reads))))
  } else {
    # Downsample full data then extract target genes
    num_downsampled_UMIs <- sapply(num_downsampled_reads, function(reads) {
      sampled_indices <- sample(reads_vec, reads)
      # Filter to target genes after downsampling
      target_data <- data[unique(sampled_indices), ] |> dplyr::filter(response_id %in% target_genes)
      return(nrow(target_data))
    })
  }

  # Return UMIs per cell
  return(num_downsampled_UMIs / cell_num)
}

# Helper function on Read-UMI table for Gasperini data
obtain_read_umi_table_Gasperini <- function(downsampling_ratio, num_SRR = 1) {
  
  qc_df <- load_gasperini_data(num_SRR)
  
  fitted_read_umi_curve <- perturbplan:::library_estimation(QC_data = qc_df)
  print("Finish the curve learning")
  
  library_size_Gasperini <- downsample_and_count_UMI(data = qc_df, ratios = downsampling_ratio)
  reads_per_cell_Gasperini <- downsampling_ratio * sum(qc_df$num_reads) / length(unique(qc_df$cell_id))
  read_umi_df <- data.frame(
    library_size = library_size_Gasperini,
    reads_per_cell = reads_per_cell_Gasperini,
    assay = ifelse(num_SRR == 1, "Gasperini (1 SRR)", sprintf("Gasperini (%d SRRs)", num_SRR))
  ) 
  print("Finish the downsampling")
  
  return(list(
    read_umi_df = read_umi_df,
    library_parameters = fitted_read_umi_curve
  ))
}

# Helper function on loading Gasperini data
load_gasperini_data <- function(num_SRR) {
  path_to_plan_data <- obtain_file_path(experiment = "Gasperini")
  SRR_dirs_list <- list.files(path = path_to_plan_data["raw"], pattern = "^SRR", full.names = TRUE)
  if (num_SRR > length(SRR_dirs_list)) {
    stop("Specified number of SRRs exceed the total number of available SRRs!")
  } else {
    SRR_dirs_in_use <- SRR_dirs_list[1:num_SRR]
  }
  
  qc_df <- NULL
  for (SRR_dir in SRR_dirs_in_use) {
    path_to_outs_folder <- sprintf("%s/outs", SRR_dir)
    raw_count_file_path <- sprintf("%s/molecule_info.h5", path_to_outs_folder)
    qc_info_file_path <- sprintf("%s/filtered_feature_bc_matrix.h5", path_to_outs_folder)
    
    raw_count_file <- rhdf5::h5read(raw_count_file_path, "count")
    umi_idx <- rhdf5::h5read(raw_count_file_path, "umi")
    
    barcode_idx <- rhdf5::h5read(raw_count_file_path, "barcode_idx")
    cell_barcodes <- rhdf5::h5read(raw_count_file_path, "barcodes")
    cell_idx <- cell_barcodes[barcode_idx + 1]
    
    gem_group <- rhdf5::h5read(raw_count_file_path, "gem_group")
    cell_id_with_gem <- paste(cell_idx, gem_group, sep = "-")
    
    RNA_idx <- rhdf5::h5read(raw_count_file_path, "feature_idx")
    gene_reference <- rhdf5::h5read(raw_count_file_path, "features")
    gene_idx <- gene_reference$id[RNA_idx + 1]
    
    raw_data_frame <- data.frame(
      num_reads = raw_count_file,
      UMI_id = umi_idx + 1,
      cell_id = cell_id_with_gem,
      response_id = gene_idx
    )
    
    qc_cell <- rhdf5::h5read(qc_info_file_path, "matrix/barcodes")
    
    qc_df <- raw_data_frame |> 
      dplyr::filter(cell_id %in% qc_cell) |>
      dplyr::mutate(cell_id = sprintf("%s_%s", cell_id, SRR_dir)) |>
      dplyr::bind_rows(qc_df)
    
    rm(qc_cell)
    print(basename(SRR_dir))
  }
  
  return(qc_df)
}

# Helper function on obtaining saturation curves varying pilot sequencing depth
obtain_read_umi_table_Gasperini_varying_sd <- function(downsampling_ratio_grid, downsampling_ratio_list, num_SRR = 1) {

  qc_df <- load_gasperini_data(num_SRR)

  # Create downsampled datasets for different sequencing depths
  qc_df_list <- vector("list", length(downsampling_ratio_grid))

  for (i in seq_along(downsampling_ratio_grid)) {
    qc_df_list[[i]] <- qc_df |>
      dplyr::mutate(
        num_reads = stats::rbinom(dplyr::n(), size = num_reads, prob = downsampling_ratio_grid[i])
      ) |>
      dplyr::filter(num_reads > 0)
  }

  library_size_Gasperini <- downsample_and_count_UMI(data = qc_df, ratios = downsampling_ratio_grid)
  reads_per_cell_Gasperini <- downsampling_ratio_grid * sum(qc_df$num_reads) / length(unique(qc_df$cell_id))
  read_umi_df <- data.frame(
    downsampling_ratio = downsampling_ratio_grid,
    library_size = library_size_Gasperini,
    reads_per_cell = reads_per_cell_Gasperini,
    assay = ifelse(num_SRR == 1, "Gasperini (1 SRR)", sprintf("Gasperini (%d SRRs)", num_SRR))
  ) |>
    dplyr::mutate(sequencing_saturation = 1 - library_size / reads_per_cell)
  print("Finish the downsampling")

  return(list(
    read_umi_df = read_umi_df,
    qc_df_list = qc_df_list
  ))
}

# Helper function on loading 10x K562 data
load_10x_data <- function() {
  # output directory
  path_to_raw_data <- .get_config_path("LOCAL_10x_K562_DATA_DIR")
  output_dir <- paste0(path_to_raw_data, "/SRRK562/outs/")

  ############################# Preprocess h5 data ###############################
  raw_path <- paste0(output_dir, list.files(output_dir, "molecule_info.h5"))
  qc_path  <- paste0(output_dir, list.files(output_dir, "feature_bc_matrix.h5"))

  # subset necessary information
  count   <- rhdf5::h5read(raw_path, "count")
  umi_idx <- rhdf5::h5read(raw_path, "umi")
  barcodes <- rhdf5::h5read(raw_path, "barcodes")
  barcode_idx <- rhdf5::h5read(raw_path, "barcode_idx")
  gem_group <- rhdf5::h5read(raw_path, "gem_group")
  cell_id <- paste(barcodes[barcode_idx + 1], gem_group, sep = "-")
  RNA_idx <- rhdf5::h5read(raw_path, "feature_idx")
  gene_id <- rhdf5::h5read(raw_path, "features")$id[RNA_idx + 1]

  # construct raw read-UMI table
  raw_df <- data.frame(num_reads = count,
                       UMI_id    = umi_idx + 1,
                       cell_id   = cell_id,
                       response_id = gene_id)

  # obtain read_umi_table
  qc_cells <- rhdf5::h5read(qc_path, "matrix/barcodes")

  # obtain the QC-ed dataframe
  qc_df <- dplyr::filter(raw_df, cell_id %in% qc_cells)

  return(qc_df)
}

# Helper function on loading T CD8 data
load_tcd8_data <- function() {
  # output directory
  path_to_raw_data <- paste0(.get_config_path("LOCAL_SHIFRUT_2018_DATA_DIR"),"/processed/")
  output_dir <- paste0(list.dirs(path_to_raw_data, recursive = FALSE, full.names = TRUE), "/outs/")

  ############################# Preprocess h5 data ###############################
  raw_path <- paste0(output_dir, list.files(output_dir, "molecule_info.h5"))
  qc_path  <- paste0(output_dir, list.files(output_dir, "filtered_feature_bc_matrix.h5"))

  # subset necessary information
  count   <- rhdf5::h5read(raw_path, "count")
  umi_idx <- rhdf5::h5read(raw_path, "umi")
  barcodes <- rhdf5::h5read(raw_path, "barcodes")
  barcode_idx <- rhdf5::h5read(raw_path, "barcode_idx")
  gem_group <- rhdf5::h5read(raw_path, "gem_group")
  cell_id <- paste(barcodes[barcode_idx + 1], gem_group, sep = "-")
  RNA_idx <- rhdf5::h5read(raw_path, "feature_idx")
  gene_id <- rhdf5::h5read(raw_path, "features")$id[RNA_idx + 1]

  # construct raw read-UMI table
  raw_df <- data.frame(num_reads = count,
                       UMI_id    = umi_idx + 1,
                       cell_id   = cell_id,
                       response_id = gene_id)

  # obtain read_umi_table
  qc_cells <- rhdf5::h5read(qc_path, "matrix/barcodes")

  # obtain the QC-ed dataframe
  qc_df <- dplyr::filter(raw_df, cell_id %in% qc_cells)

  return(qc_df)
}

# Helper function on Read-UMI table for T CD8 cell type
obtain_read_umi_table_tcd8 <- function(downsampling_ratio) {

  qc_df <- load_tcd8_data()

  # downsample to get library size
  library_size_tcd8 <- downsample_and_count_UMI(data = qc_df, ratios = downsampling_ratio)
  reads_per_cell_tcd8 <- downsampling_ratio * sum(qc_df$num_reads) / length(unique(qc_df$cell_id))
  read_umi_df <- data.frame(
    library_size = library_size_tcd8,
    reads_per_cell = reads_per_cell_tcd8,
    assay = "TCD8 example"
  )

  # return the QC data
  return(read_umi_df)
}

# Helper function on obtaining saturation curves varying pilot sequencing depth for T CD8 data
obtain_read_umi_table_tcd8_varying_sd <- function(downsampling_ratio_grid, downsampling_ratio_list) {

  qc_df <- load_tcd8_data()

  # Create downsampled datasets for different sequencing depths
  qc_df_list <- vector("list", length(downsampling_ratio_grid))

  for (i in seq_along(downsampling_ratio_grid)) {
    qc_df_list[[i]] <- qc_df |>
      dplyr::mutate(
        num_reads = stats::rbinom(dplyr::n(), size = num_reads, prob = downsampling_ratio_grid[i])
      ) |>
      dplyr::filter(num_reads > 0)
  }

  library_size_tcd8 <- downsample_and_count_UMI(data = qc_df, ratios = downsampling_ratio_grid)
  reads_per_cell_tcd8 <- downsampling_ratio_grid * sum(qc_df$num_reads) / length(unique(qc_df$cell_id))
  read_umi_df <- data.frame(
    downsampling_ratio = downsampling_ratio_grid,
    library_size = library_size_tcd8,
    reads_per_cell = reads_per_cell_tcd8,
    assay = "TCD8"
  ) |>
    dplyr::mutate(sequencing_saturation = 1 - library_size / reads_per_cell)
  print("Finish the downsampling")

  return(list(
    read_umi_df = read_umi_df,
    qc_df_list = qc_df_list
  ))
}

# Helper function on obtaining Read-UMI table for 10x data
obtain_read_umi_table_10x <- function(downsampling_ratio) {

  qc_df <- load_10x_data()

  # downsample to get library size
  library_size_10x <- downsample_and_count_UMI(data = qc_df, ratios = downsampling_ratio)
  reads_per_cell_10x <- downsampling_ratio * sum(qc_df$num_reads) / length(unique(qc_df$cell_id))
  read_umi_df <- data.frame(
    library_size = library_size_10x,
    reads_per_cell = reads_per_cell_10x,
    assay = "10x example"
  )

  # return the QC data
  return(read_umi_df)
}

# Helper function on obtaining saturation curves varying pilot sequencing depth for 10x data
obtain_read_umi_table_10x_varying_sd <- function(downsampling_ratio_grid, downsampling_ratio_list) {

  qc_df <- load_10x_data()

  # Create downsampled datasets for different sequencing depths
  qc_df_list <- vector("list", length(downsampling_ratio_grid))

  for (i in seq_along(downsampling_ratio_grid)) {
    qc_df_list[[i]] <- qc_df |>
      dplyr::mutate(
        num_reads = stats::rbinom(dplyr::n(), size = num_reads, prob = downsampling_ratio_grid[i])
      ) |>
      dplyr::filter(num_reads > 0)
  }

  library_size_10x <- downsample_and_count_UMI(data = qc_df, ratios = downsampling_ratio_grid)
  reads_per_cell_10x <- downsampling_ratio_grid * sum(qc_df$num_reads) / length(unique(qc_df$cell_id))
  read_umi_df <- data.frame(
    downsampling_ratio = downsampling_ratio_grid,
    library_size = library_size_10x,
    reads_per_cell = reads_per_cell_10x,
    assay = "10x K562"
  ) |>
    dplyr::mutate(sequencing_saturation = 1 - library_size / reads_per_cell)
  print("Finish the downsampling")

  return(list(
    read_umi_df = read_umi_df,
    qc_df_list = qc_df_list
  ))
}

# Helper function on obtaining saturation curves varying pilot sequencing depth for K562 TAP data
obtain_read_umi_table_k562_tap_varying_sd <- function(downsampling_ratio_grid, downsampling_ratio_list) {

  # Use K562_Ray data from perturbplan package
  data("K562_Ray", package = "perturbplan")

  # output directory from config
  path_to_raw_data <- .get_config_path("LOCAL_RAY_2025_RAW_DATA_DIR")
  output_dir <- paste0(path_to_raw_data, "/processed/perturbplan-demo/outs/")

  ######################### Load target gene list ################################
  # Read target_panel.csv to get the gene list
  panel_path <- file.path(output_dir, "target_panel.csv")
  if (!file.exists(panel_path)) {
    stop("target_panel.csv not found at: ", panel_path)
  }

  # helper: safe trim
  .trim <- function(x) gsub("^\\s+|\\s+$", "", x)

  gene_list <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    # first column, no header, start from row 7 -> skip first 6 lines
    dt <- data.table::fread(panel_path, header = FALSE, skip = 6, select = 1L, data.table = FALSE)
    gene_list <- unique(.trim(dt[[1]]))
  } else {
    # base R fallback
    con <- file(panel_path, open = "r"); on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    if (length(lines) < 7) stop("target_panel.csv has fewer than 7 lines.")
    first_col <- vapply(strsplit(lines[-(1:6)], ",", fixed = TRUE), function(z) z[[1]], character(1L))
    gene_list <- unique(.trim(first_col))
  }

  # drop empties / NAs just in case
  gene_list <- gene_list[!is.na(gene_list) & nzchar(gene_list)]

  if (length(gene_list) == 0L) {
    stop("Parsed gene_list is empty from ", panel_path)
  }

  ############################# Preprocess h5 data ###############################
  raw_path <- paste0(output_dir, list.files(output_dir, "molecule_info.h5"))
  qc_path  <- paste0(output_dir, list.files(output_dir, "filtered_feature_bc_matrix.h5"))

  # subset necessary information
  count   <- rhdf5::h5read(raw_path, "count")
  umi_idx <- rhdf5::h5read(raw_path, "umi")
  barcodes <- rhdf5::h5read(raw_path, "barcodes")
  barcode_idx <- rhdf5::h5read(raw_path, "barcode_idx")
  gem_group <- rhdf5::h5read(raw_path, "gem_group")
  cell_id <- paste(barcodes[barcode_idx + 1], gem_group, sep = "-")
  RNA_idx <- rhdf5::h5read(raw_path, "feature_idx")
  gene_id <- rhdf5::h5read(raw_path, "features")$id[RNA_idx + 1]

  # construct raw read-UMI table
  raw_df <- data.frame(num_reads = count,
                       UMI_id    = umi_idx + 1,
                       cell_id   = cell_id,
                       response_id = gene_id)

  # obtain read_umi_table
  qc_cells <- rhdf5::h5read(qc_path, "matrix/barcodes")

  # obtain the QC-ed dataframe and filter by target gene list
  qc_df <- raw_df |>
    dplyr::filter(cell_id %in% qc_cells) |>
    dplyr::filter(response_id %in% gene_list)

  # Create downsampled datasets for different sequencing depths
  qc_df_list <- vector("list", length(downsampling_ratio_grid))

  for (i in seq_along(downsampling_ratio_grid)) {
    qc_df_list[[i]] <- qc_df |>
      dplyr::mutate(
        num_reads = stats::rbinom(dplyr::n(), size = num_reads, prob = downsampling_ratio_grid[i])
      ) |>
      dplyr::filter(num_reads > 0)
  }

  library_size_tap <- downsample_and_count_UMI(data = qc_df, ratios = downsampling_ratio_grid)
  reads_per_cell_tap <- downsampling_ratio_grid * sum(qc_df$num_reads) / length(unique(qc_df$cell_id))
  read_umi_df <- data.frame(
    downsampling_ratio = downsampling_ratio_grid,
    library_size = library_size_tap,
    reads_per_cell = reads_per_cell_tap,
    assay = "K562 TAP-seq"
  ) |>
    dplyr::mutate(sequencing_saturation = 1 - library_size / reads_per_cell)
  print("Finish the downsampling for K562 TAP-seq")

  return(list(
    read_umi_df = read_umi_df,
    qc_df_list = qc_df_list
  ))
}

# -----------------------------
# Alternative library estimation method using downsampling
# -----------------------------
library_estimation_alternative <- function(QC_data,
                                          downsample_ratio = c(0.3, 0.5, 0.7),
                                          D2_rough = c(0.3, 0.5, 0.7)) {
  result <- library_computation_alternative(QC_data, downsample_ratio, D2_rough)

  # Check if fitted_model exists and is valid
  if (is.null(result$fitted_model)) {
    stop("library_computation_alternative returned NULL fitted_model")
  }

  # Extract coefficients
  model_coefs <- stats::coef(result$fitted_model)

  # Check if coefficients exist
  if (is.null(model_coefs) || length(model_coefs) == 0) {
    stop("Failed to extract coefficients from fitted model")
  }

  total_UMIs <- model_coefs["total_UMIs"]
  umi_variation <- model_coefs["D2"]

  # Check if coefficients were found
  if (is.na(total_UMIs)) {
    stop("total_UMIs coefficient not found in model. Available coefficients: ",
         paste(names(model_coefs), collapse = ", "))
  }
  if (is.na(umi_variation)) {
    stop("D2 coefficient not found in model. Available coefficients: ",
         paste(names(model_coefs), collapse = ", "))
  }

  return(list(
    UMI_per_cell = unname(as.numeric(total_UMIs)),
    variation = unname(as.numeric(umi_variation)),
    best_downsample_ratio = result$best_downsample_ratio,
    best_D2_rough = result$best_D2_rough,
    best_relative_error = result$best_relative_error
  ))
}

library_computation_alternative <- function(QC_data,
                                            downsample_ratio = c(0.3, 0.5, 0.7),
                                            D2_rough = c(0.3, 0.5, 0.7),
                                            seed = 42) {

  # Pre-compute data that's common across all combinations
  cell_num <- length(unique(QC_data$cell_id))
  num_observed_reads <- sum(QC_data$num_reads)
  reads_vec <- rep(seq_along(QC_data$num_reads), QC_data$num_reads)
  num_observed_umis <- length(unique(reads_vec))

  # Create all combinations of downsample_ratio and D2_rough
  param_grid <- expand.grid(
    downsample_ratio = downsample_ratio,
    D2_rough = D2_rough,
    stringsAsFactors = FALSE
  )

  # Store all results
  all_results <- list()
  result_counter <- 1

  # Try each combination of downsample_ratio and D2_rough
  for (i in seq_len(nrow(param_grid))) {
    ds_ratio <- param_grid$downsample_ratio[i]
    d2_rough <- param_grid$D2_rough[i]

    ########################### downsample the data ##############################
    # Set seed for reproducibility (different seed for each parameter combination)
    if (!is.null(seed)) {
      set.seed(seed + i)
    }

    # perform downsampling and append the results together with observed reads-UMIs
    num_downsampled_reads <- round(num_observed_reads * ds_ratio)
    num_downsampled_UMIs <- length(unique(sample(reads_vec, num_downsampled_reads)))
    down_sample_added <- data.frame(
      num_reads = num_downsampled_reads / cell_num,
      num_UMIs = num_downsampled_UMIs / cell_num
    )
    down_sample_df <- down_sample_added |>
      dplyr::bind_rows(data.frame(
        num_reads = num_observed_reads / cell_num,
        num_UMIs = num_observed_umis / cell_num
      )) |>
      dplyr::arrange(num_reads, num_UMIs)

    ####################### fit nonlinear model ##################################
    # Calculate two initial values for total_UMIs
    delicate_initial <- (1 + d2_rough) * (num_observed_reads / cell_num)^2 /
      (2 * (num_observed_reads - num_observed_umis) / cell_num)
    rough_initial <- num_observed_umis / cell_num
    inital_num_UMIs_vec <- stats::setNames(
      c(delicate_initial, rough_initial),
      c("delicate", "rough")
    )

    # fit model with different initial values on total UMIs
    for (init_name in names(inital_num_UMIs_vec)) {
      initial_UMIs <- unname(inital_num_UMIs_vec[init_name])

      # do the model fitting with error handling
      nlm_fitting <- tryCatch({
        minpack.lm::nlsLM(
          num_UMIs ~ total_UMIs * (1 - exp(-num_reads / total_UMIs) *
                                      (1 + D2 * num_reads^2 / (2 * total_UMIs^2))),
          data = down_sample_df,
          start = list(total_UMIs = initial_UMIs, D2 = d2_rough),
          upper = c(Inf, 1),
          lower = c(0, 0),
          control = minpack.lm::nls.lm.control(maxiter = 200, ftol = 1e-6, ptol = 1e-6)
        )
      }, error = function(e) {
        return(NULL)
      })

      if (!is.null(nlm_fitting)) {
        # Calculate relative loss
        relative_loss <- sum((stats::predict(nlm_fitting) / (down_sample_df$num_UMIs) - 1)^2)

        # Store result
        all_results[[result_counter]] <- list(
          fitted_model = nlm_fitting,
          relative_error = relative_loss,
          downsample_ratio = ds_ratio,
          D2_rough = d2_rough,
          init_method = init_name
        )
        result_counter <- result_counter + 1
      }
    }
  }

  # Check if we have any valid fits
  if (length(all_results) == 0) {
    stop("All nlsLM fitting attempts failed across all parameter combinations. ",
         "Consider using a different method or checking your data.")
  }

  # Select the best fit (minimum relative error)
  relative_errors <- sapply(all_results, function(x) x$relative_error)
  best_idx <- which.min(relative_errors)
  best_result <- all_results[[best_idx]]

  # Add a warning about the relative error
  if (!is.null(best_result$relative_error) &&
      !is.na(best_result$relative_error) &&
      best_result$relative_error > 0.05) {
    perc_error <- round(100 * best_result$relative_error, 2)
    message(
      sprintf("Best fit relative error: %.2f%% (downsample_ratio=%.1f, D2_rough=%.1f, init=%s)",
              perc_error,
              best_result$downsample_ratio,
              best_result$D2_rough,
              best_result$init_method)
    )
  }

  return(list(
    fitted_model = best_result$fitted_model,
    best_downsample_ratio = best_result$downsample_ratio,
    best_D2_rough = best_result$D2_rough,
    best_relative_error = best_result$relative_error,
    best_init_method = best_result$init_method
  ))
}

# Prediction function for alternative method
fit_read_UMI_curve_alternative <- function(reads_per_cell, UMI_per_cell, variation) {
  UMI_per_cell * (1 - exp(-reads_per_cell / UMI_per_cell) * (1 + variation * reads_per_cell^2 / (2 * UMI_per_cell^2)))
}

library_estimation_preseqR_f <- function(QC_data, mt = 20) {
  # Filter out NA values before creating histogram
  valid_reads <- QC_data$num_reads[!is.na(QC_data$num_reads) & QC_data$num_reads > 0]

  if (length(valid_reads) == 0) {
    stop("No valid read counts found in QC_data after filtering NAs and zeros")
  }

  read_umi_summary <- table(valid_reads)
  num_cells <- length(unique(QC_data$cell_id))

  # Extract read counts (names) and frequencies (values)
  preseq_input <- cbind(as.integer(names(read_umi_summary)), as.vector(read_umi_summary))

  # Additional validation: remove any rows with NA
  preseq_input <- preseq_input[complete.cases(preseq_input), , drop = FALSE]

  if (nrow(preseq_input) == 0) {
    stop("preseq_input has no valid rows after filtering")
  }

  reads_per_cell_original <- sum(preseq_input[, 1] * preseq_input[, 2]) / num_cells

  # Follow preseqR.rSAC logic exactly
  para <- preseqR::preseqR.ztnb.em(preseq_input)
  shape <- para$size
  mu <- para$mu

  if (shape <= 1) {
    # Use ds.rSAC method (RFA)
    method_used <- "RFA"

    # Call ds.rSAC following preseqR.rSAC logic
    preseqR_fn <- suppressWarnings(preseqR::ds.rSAC(n = preseq_input, mt = mt))
    fn_env <- environment(preseqR_fn)

    # Extract parameters from ds.rSAC environment
    coefs <- fn_env$coefs
    poles <- fn_env$poles
    valid_estimator <- fn_env$valid.estimator

    if (!is.null(valid_estimator) && valid_estimator == FALSE) {
      # Invalid estimator - will return constant
      n_data <- fn_env$n
      UMI_per_cell_at_saturation <- sum(n_data[, 2]) / num_cells

      # Store parameters for invalid estimator case
      params <- list(
        method_used = method_used,
        valid_estimator = FALSE,
        constant_value = as.numeric(sum(n_data[, 2])),
        reads_norm = as.numeric(reads_per_cell_original),
        n_cells = as.numeric(num_cells),
        UMI_per_cell_at_saturation = UMI_per_cell_at_saturation
      )
    } else {
      # Valid RFA estimator
      UMI_per_cell_at_saturation <- as.numeric(Re(sum(coefs))) / num_cells

      # Store complex numbers as real and imaginary parts
      params <- list(
        method_used = method_used,
        valid_estimator = TRUE,
        coefs_real = as.numeric(Re(coefs)),
        coefs_imag = as.numeric(Im(coefs)),
        poles_real = as.numeric(Re(poles)),
        poles_imag = as.numeric(Im(poles)),
        reads_norm = as.numeric(reads_per_cell_original),
        n_cells = as.numeric(num_cells),
        UMI_per_cell_at_saturation = UMI_per_cell_at_saturation
      )
    }
  } else {
    # Use ZTNB closed-form method
    method_used <- "ZTNB"

    # Follow the exact formula from preseqR.rSAC
    p <- 1 - dnbinom(0, size = shape, mu = mu)
    L <- sum(as.numeric(preseq_input[, 2])) / p

    # At saturation (infinite reads), the ZTNB formula approaches L
    UMI_per_cell_at_saturation <- L / num_cells

    # Store parameters needed for ZTNB formula
    params <- list(
      method_used = method_used,
      L = as.numeric(L),
      size = as.numeric(shape),
      mu = as.numeric(mu),
      reads_norm = as.numeric(reads_per_cell_original),
      n_cells = as.numeric(num_cells),
      UMI_per_cell_at_saturation = UMI_per_cell_at_saturation
    )
  }

  params
}

fit_read_UMI_curve_preseqR_f <- function(reads_per_cell, rSAC_fn_wrapper){
  # Validate inputs
  if (length(reads_per_cell) == 0) {
    stop("reads_per_cell has length 0")
  }

  # Extract normalization parameters and method
  reads_norm <- rSAC_fn_wrapper$reads_norm
  n_cells <- rSAC_fn_wrapper$n_cells
  method_used <- rSAC_fn_wrapper$method_used

  # Normalize reads_per_cell to the scale used during estimation
  t <- reads_per_cell / reads_norm

  if (method_used == "ZTNB") {
    # Use ZTNB closed-form formula: L * pnbinom(r - 1, size = size, mu = mu * t, lower.tail = FALSE)
    L <- rSAC_fn_wrapper$L
    size <- rSAC_fn_wrapper$size
    mu <- rSAC_fn_wrapper$mu

    predictions <- sapply(t, function(x) {
      L * pnbinom(0, size = size, mu = mu * x, lower.tail = FALSE)
    }) / n_cells

  } else if (method_used == "RFA") {
    # Use RFA (ds.rSAC) formula
    if (!rSAC_fn_wrapper$valid_estimator) {
      # Invalid estimator - return constant: sum(n[, 2])
      predictions <- rep(rSAC_fn_wrapper$constant_value / n_cells, length(t))
    } else {
      # Valid RFA estimator - use formula: Re(coefs %*% (x/(x - poles))^r)
      # Reconstruct complex coefficients and poles
      coefs <- complex(real = rSAC_fn_wrapper$coefs_real,
                      imaginary = rSAC_fn_wrapper$coefs_imag)
      poles <- complex(real = rSAC_fn_wrapper$poles_real,
                      imaginary = rSAC_fn_wrapper$poles_imag)

      # Apply the ds.rSAC formula
      predictions <- sapply(t, function(x) {
        Re(coefs %*% (x/(x - poles)))
      }) / n_cells
    }
  } else {
    stop(sprintf("Unknown method_used: %s", method_used))
  }

  # Handle edge cases: ensure no negative values, NaN, or Inf
  predictions[is.na(predictions)] <- 0
  predictions[is.infinite(predictions)] <- 0
  predictions[predictions < 0] <- 0

  return(predictions)
}

# Pure ds.rSAC method (without ZTNB fallback)
library_estimation_ds_rSAC <- function(QC_data, r = 1, mt = 20) {
  # Filter out NA values before creating histogram
  valid_reads <- QC_data$num_reads[!is.na(QC_data$num_reads) & QC_data$num_reads > 0]

  if (length(valid_reads) == 0) {
    stop("No valid read counts found in QC_data after filtering NAs and zeros")
  }

  read_umi_summary <- table(valid_reads)
  num_cells <- length(unique(QC_data$cell_id))

  # Extract read counts (names) and frequencies (values)
  preseq_input <- cbind(as.integer(names(read_umi_summary)), as.vector(read_umi_summary))

  # Additional validation: remove any rows with NA
  preseq_input <- preseq_input[complete.cases(preseq_input), , drop = FALSE]

  if (nrow(preseq_input) == 0) {
    stop("preseq_input has no valid rows after filtering")
  }

  reads_per_cell_original <- sum(preseq_input[, 1] * preseq_input[, 2]) / num_cells

  # Directly call ds.rSAC without ZTNB
  # Suppress polynomial warnings from preseqR internal conversions
  preseqR_fn <- suppressWarnings(preseqR::ds.rSAC(n = preseq_input, r = r, mt = mt))

  # Extract UMI_per_cell from preseqR function environment
  fn_env <- environment(preseqR_fn)
  coefs <- fn_env$coefs
  valid_estimator <- fn_env$valid.estimator

  # Check if coefficients have significant imaginary parts
  if (!is.null(coefs) && is.complex(coefs)) {
    max_imag <- max(abs(Im(coefs)))
    max_real <- max(abs(Re(coefs)))
    if (max_imag > 1e-10 && max_imag / max_real > 1e-6) {
      message(sprintf("Note: ds.rSAC coefficients have imaginary parts (max |Im|/|Re| ratio: %.2e)",
                      max_imag / max_real))
    }
  }

  # Check if this is a valid estimator
  if (!is.null(valid_estimator) && valid_estimator == FALSE) {
    # Invalid estimator
    warning("ds.rSAC returned invalid estimator (insufficient data depth). Using constant prediction.")
    n_data <- fn_env$n
    UMI_per_cell_at_saturation <- sum(n_data[, 2]) / num_cells
  } else if (is.null(coefs)) {
    # No coefficients found
    stop("ds.rSAC returned function without coefficients")
  } else {
    # Valid ds.rSAC estimator
    UMI_per_cell_at_saturation <- as.numeric(Re(sum(coefs))) / num_cells
  }

  # Clean up the function environment to reduce file size
  # Remove large intermediate objects while keeping necessary parameters
  fn_env <- environment(preseqR_fn)

  # List all objects in the environment
  all_vars <- ls(envir = fn_env, all.names = TRUE)

  # Define essential variables to keep for ds.rSAC method
  # Note: "n" is NOT needed - we've already computed UMI_per_cell_at_saturation
  essential_vars <- c("coefs", "poles", "r", "valid.estimator")

  # Remove non-essential variables
  vars_to_remove <- setdiff(all_vars, essential_vars)
  if (length(vars_to_remove) > 0) {
    rm(list = vars_to_remove, envir = fn_env)
  }

  # Store the cleaned function with normalization metadata
  preseqR_fn_wrapper <- list(
    fn = preseqR_fn,
    reads_norm = as.numeric(reads_per_cell_original),
    n_cells = as.numeric(num_cells)
  )

  # Note: ds.rSAC doesn't provide size parameter, so we can't compute variation
  # Return function wrapper instead of manually extracted parameters
  list(
    rSAC_fn = preseqR_fn_wrapper,
    UMI_per_cell_at_saturation = UMI_per_cell_at_saturation,
    UMI_richness_variation = NA  # Not available for pure ds.rSAC
  )
}

fit_read_UMI_curve_ds_rSAC <- function(reads_per_cell, rSAC_fn_wrapper){
  # Validate inputs
  if (length(reads_per_cell) == 0) {
    stop("reads_per_cell has length 0")
  }

  # Extract function and normalization parameters
  preseqR_fn <- rSAC_fn_wrapper$fn
  reads_norm <- rSAC_fn_wrapper$reads_norm
  n_cells <- rSAC_fn_wrapper$n_cells

  # Normalize reads_per_cell to the scale used during estimation
  t <- reads_per_cell / reads_norm

  # Use preseqR's own function for prediction
  predictions <- preseqR_fn(t) / n_cells

  # Handle edge cases: ensure no negative values, NaN, or Inf
  predictions[is.na(predictions)] <- 0
  predictions[is.infinite(predictions)] <- 0
  predictions[predictions < 0] <- 0

  return(predictions)
}



# -----------------------------
# IO helpers (one output, elegant)
# -----------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_plot <- function(plot, path, width, height) {
  ensure_dir(dirname(path))
  ggsave(filename = path, plot = plot, width = width, height = height)
  invisible(path)
}