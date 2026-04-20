library(perturbplan)
library(dplyr)
library(tidyr)
library(tibble)
library(minpack.lm)
library(preseqR)

# set seed and source necessary files
set.seed(20)
source("~/.Rprofile")
source("reproduction-code-prep/saturation-curve-fitting/helper.R")

# specify the directory to results
dir_to_results <- "reproduction-code-prep/saturation-curve-fitting/results"

# -----------------------------
# Global configs
# -----------------------------
downsampling_ratio_list <- 10^{seq(-2, 0, length.out = 21)}
reads_per_cell_grid_validation <- 10^{seq(2, 4.5, length.out = 500)}
reads_per_cell_grid_varying_sd      <- 10^{seq(2, 6, length.out = 500)}
downsampling_ratio_grid_varying_sd  <- c(0.01, 0.03, 0.1, 0.3, 1)
num_SRR_in_use <- 1

# -----------------------------
# Utilities
# -----------------------------
load_or_compute <- function(read_cache, param_cache, compute_fn, post_read_fn = NULL,
                            qc_cache = NULL, qc_df_list = NULL) {
  # If qc_df_list is provided, use it directly (shared downsampling)
  if (!is.null(qc_df_list)) {
    if (file.exists(param_cache)) {
      params <- readRDS(param_cache)
      read_df <- readRDS(read_cache)
    } else {
      res <- compute_fn(qc_df_list)
      read_df <- res$read_umi_df
      params  <- res$library_parameters

      # cache
      saveRDS(read_df, read_cache)
      saveRDS(params,  param_cache)
    }
    if (!is.null(post_read_fn)) read_df <- post_read_fn(read_df)
    return(list(read_umi_df = read_df, library_parameters = params))
  }

  # Original behavior: load from cache or compute
  if (file.exists(read_cache) && file.exists(param_cache)) {
    read_df <- readRDS(read_cache)
    params  <- readRDS(param_cache)
    # Load qc_df_list if requested
    if (!is.null(qc_cache) && file.exists(qc_cache)) {
      qc_list <- readRDS(qc_cache)
    } else {
      qc_list <- NULL
    }
  } else {
    res <- compute_fn()
    read_df <- res$read_umi_df
    params  <- res$library_parameters
    qc_list <- res$qc_df_list

    # cache (speed-up; NOT an output artifact)
    saveRDS(read_df, read_cache)
    saveRDS(params,  param_cache)
    if (!is.null(qc_cache)) {
      saveRDS(qc_list, qc_cache)
    }
  }
  if (!is.null(post_read_fn)) read_df <- post_read_fn(read_df)
  list(read_umi_df = read_df, library_parameters = params, qc_df_list = qc_list)
}

fit_log_model_quantiles <- function(df, probs = c(0.25, 0.50, 0.75, 1.00)) {
  n <- nrow(df)
  idx <- unique(pmax(1, pmin(n, round(n * probs))))
  lm(library_size ~ log10(reads_per_cell), data = df[idx, , drop = FALSE])
}

make_validation_grid <- function(read_umi_df, library_parameters, log_model,
                                 cell_type, reads_per_cell_grid) {
  tibble(
    reads_per_cell = reads_per_cell_grid,
    library_size   = perturbplan:::fit_read_UMI_curve_cpp(
      reads_per_cell = reads_per_cell_grid,
      rSAC_fn_wrapper = library_parameters
    ),
    source    = "Perturbplan fit",
    cell_type = cell_type
  ) %>%
    bind_rows(
      read_umi_df %>%
        transmute(
          reads_per_cell, library_size,
          source    = "Downsampled real data",
          cell_type = cell_type
        )
    ) %>%
    bind_rows(
      tibble(
        reads_per_cell = reads_per_cell_grid,
        library_size   = predict(log_model, newdata = data.frame(reads_per_cell = reads_per_cell_grid)),
        source         = "scPower fit",
        cell_type      = cell_type
      )
    )
}

# ============================================================
# Part A: Validation plots (T_CD8 + K562)  -> panels a/b
# ============================================================

# ---- T_CD8_Shifrut ----
data("T_CD8_Shifrut")

# construct Read-UMI table
read_umi_path_tcd8 <- "reproduction-code-prep/saturation-curve-fitting/results/read_umi_df_tcd8.rds"
if (file.exists(read_umi_path_tcd8)) {
  read_umi_df_tcd8 <- readRDS(read_umi_path_tcd8)
} else {
  read_umi_df_tcd8 <- obtain_read_umi_table_tcd8(downsampling_ratio = downsampling_ratio_list)
  saveRDS(read_umi_df_tcd8, read_umi_path_tcd8) # cache only
}

# fit log-linear model
log_model_tcd8 <- fit_log_model_quantiles(read_umi_df_tcd8)

# fit library parameters for T_CD8
library_parameters_tcd8 <- T_CD8_Shifrut$library_parameters

# make validation grid
grid_tcd8 <- make_validation_grid(
  read_umi_df         = read_umi_df_tcd8,
  library_parameters  = library_parameters_tcd8,
  log_model           = log_model_tcd8,
  cell_type           = "T CD8 (Shifrut)",
  reads_per_cell_grid = reads_per_cell_grid_validation
)

# ---- K562_Gasperini (standard) ----
data("K562_Gasperini")
read_umi_cache_k562 <- "reproduction-code-prep/saturation-curve-fitting/results/read_umi_df_k562.rds"
param_cache_k562    <- sprintf("reproduction-code-prep/saturation-curve-fitting/results/library_parameters_on_%d_SRRs.rds", num_SRR_in_use)

# Construct Read-UMI table using fixed downsampling ratios (same as T_CD8)
k562_result <- load_or_compute(
  read_cache  = read_umi_cache_k562,
  param_cache = param_cache_k562,
  compute_fn  = function() {
    obtain_read_umi_table_Gasperini(
      downsampling_ratio = downsampling_ratio_list,
      num_SRR = num_SRR_in_use
    )
  }
)
read_umi_df_k562 <- k562_result$read_umi_df
library_parameters_k562 <- k562_result$library_parameters

# Fit log-linear model
log_model_k562 <- fit_log_model_quantiles(read_umi_df_k562)

# make validation grid
grid_k562 <- make_validation_grid(
  read_umi_df         = read_umi_df_k562,
  library_parameters  = library_parameters_k562,
  log_model           = log_model_k562,
  cell_type           = "K562 (Gasperini)",
  reads_per_cell_grid = reads_per_cell_grid_validation
)

# combine two data frames
plot_data_validation <- bind_rows(grid_tcd8, grid_k562) %>% dplyr::filter(library_size > 0)

# save data to compare two methods
saveRDS(plot_data_validation, sprintf("%s/saturation_curve_comparison.rds", dir_to_results))

# ============================================================
# Part B: Varying sequencing saturation plot (K562 varying_sd) -> panel c
# ============================================================
read_umi_cache_vsd <- "reproduction-code-prep/saturation-curve-fitting/results/read_umi_df_varying_sd.rds"
param_cache_vsd    <- sprintf("reproduction-code-prep/saturation-curve-fitting/results/library_parameters_on_%d_SRRs_varying_sd.rds", num_SRR_in_use)

# Construct Read-UMI table by varying sequencing depth
k562_vsd <- load_or_compute(
  read_cache  = read_umi_cache_vsd,
  param_cache = param_cache_vsd,
  compute_fn  = function() {
    result <- obtain_read_umi_table_Gasperini_varying_sd(
      downsampling_ratio_grid = downsampling_ratio_grid_varying_sd,
      downsampling_ratio_list = downsampling_ratio_list,
      num_SRR = num_SRR_in_use
    )

    # Use standard library_estimation method
    library_parameters <- lapply(result$qc_df_list, function(qc_df) {
      perturbplan:::library_estimation(QC_data = qc_df)
    })

    list(
      read_umi_df = result$read_umi_df,
      library_parameters = library_parameters
    )
  }
)
read_umi_df_vsd         <- k562_vsd$read_umi_df
library_parameters_vsd  <- k562_vsd$library_parameters

# Construct Read-UMI table grid
read_umi_grid_df_vsd <- read_umi_df_vsd %>%
  transmute(sequencing_saturation, downsampling_ratio, param_idx = row_number()) %>%
  crossing(reads_per_cell = reads_per_cell_grid_varying_sd) %>%
  rowwise() %>%
  mutate(
    library_size = perturbplan:::fit_read_UMI_curve_cpp(
      reads_per_cell = reads_per_cell,
      rSAC_fn_wrapper = library_parameters_vsd[[param_idx]]
    )
  ) %>%
  ungroup()

# save Read-UMI table
saveRDS(read_umi_grid_df_vsd, sprintf("%s/saturation_curve_comparison_varying_sequencing_depth.rds", dir_to_results))

# ============================================================
# Part C: Varying sequencing saturation plot (K562_10x varying_sd) -> panel d
# ============================================================
read_umi_cache_vsd_10x <- "reproduction-code-prep/saturation-curve-fitting/results/read_umi_df_varying_sd_10x.rds"
param_cache_vsd_10x    <- "reproduction-code-prep/saturation-curve-fitting/results/library_parameters_10x_varying_sd.rds"

# Construct Read-UMI table by varying sequencing depth for 10x data
k562_10x_vsd <- load_or_compute(
  read_cache  = read_umi_cache_vsd_10x,
  param_cache = param_cache_vsd_10x,
  compute_fn  = function() {
    result <- obtain_read_umi_table_10x_varying_sd(
      downsampling_ratio_grid = downsampling_ratio_grid_varying_sd,
      downsampling_ratio_list = downsampling_ratio_list
    )

    # Use standard library_estimation method
    library_parameters <- lapply(result$qc_df_list, function(qc_df) {
      perturbplan:::library_estimation(QC_data = qc_df)
    })

    list(
      read_umi_df = result$read_umi_df,
      library_parameters = library_parameters
    )
  }
)
read_umi_df_vsd_10x         <- k562_10x_vsd$read_umi_df
library_parameters_vsd_10x  <- k562_10x_vsd$library_parameters

# Construct Read-UMI table grid
read_umi_grid_df_vsd_10x <- read_umi_df_vsd_10x %>%
  transmute(sequencing_saturation, downsampling_ratio, param_idx = row_number()) %>%
  crossing(reads_per_cell = reads_per_cell_grid_varying_sd) %>%
  rowwise() %>%
  mutate(
    library_size = perturbplan:::fit_read_UMI_curve_cpp(
      reads_per_cell = reads_per_cell,
      rSAC_fn_wrapper = library_parameters_vsd_10x[[param_idx]]
    )
  ) %>%
  ungroup()

# save Read-UMI table
saveRDS(read_umi_df_vsd_10x, sprintf("%s/read_umi_df_varying_sd_10x.rds", dir_to_results))
saveRDS(library_parameters_vsd_10x, sprintf("%s/library_parameters_10x_varying_sd.rds", dir_to_results))
saveRDS(read_umi_grid_df_vsd_10x, sprintf("%s/saturation_curve_comparison_varying_sequencing_depth_10x.rds", dir_to_results))

# ============================================================
# Part D: Varying sequencing saturation plot (TCD8 varying_sd) -> panel e
# ============================================================
read_umi_cache_vsd_tcd8 <- "reproduction-code-prep/saturation-curve-fitting/results/read_umi_df_varying_sd_tcd8.rds"
param_cache_vsd_tcd8    <- "reproduction-code-prep/saturation-curve-fitting/results/library_parameters_tcd8_varying_sd.rds"

# Construct Read-UMI table by varying sequencing depth for TCD8 data
tcd8_vsd <- load_or_compute(
  read_cache  = read_umi_cache_vsd_tcd8,
  param_cache = param_cache_vsd_tcd8,
  compute_fn  = function() {
    result <- obtain_read_umi_table_tcd8_varying_sd(
      downsampling_ratio_grid = downsampling_ratio_grid_varying_sd,
      downsampling_ratio_list = downsampling_ratio_list
    )

    # Use standard library_estimation method
    library_parameters <- lapply(result$qc_df_list, function(qc_df) {
      perturbplan:::library_estimation(QC_data = qc_df)
    })

    list(
      read_umi_df = result$read_umi_df,
      library_parameters = library_parameters
    )
  }
)
read_umi_df_vsd_tcd8         <- tcd8_vsd$read_umi_df
library_parameters_vsd_tcd8  <- tcd8_vsd$library_parameters

# Construct Read-UMI table grid
read_umi_grid_df_vsd_tcd8 <- read_umi_df_vsd_tcd8 %>%
  transmute(sequencing_saturation, downsampling_ratio, param_idx = row_number()) %>%
  crossing(reads_per_cell = reads_per_cell_grid_varying_sd) %>%
  rowwise() %>%
  mutate(
    library_size = perturbplan:::fit_read_UMI_curve_cpp(
      reads_per_cell = reads_per_cell,
      rSAC_fn_wrapper = library_parameters_vsd_tcd8[[param_idx]]
    )
  ) %>%
  ungroup()

# save Read-UMI table
saveRDS(read_umi_df_vsd_tcd8, sprintf("%s/read_umi_df_varying_sd_tcd8.rds", dir_to_results))
saveRDS(library_parameters_vsd_tcd8, sprintf("%s/library_parameters_tcd8_varying_sd.rds", dir_to_results))
saveRDS(read_umi_grid_df_vsd_tcd8, sprintf("%s/saturation_curve_comparison_varying_sequencing_depth_tcd8.rds", dir_to_results))

