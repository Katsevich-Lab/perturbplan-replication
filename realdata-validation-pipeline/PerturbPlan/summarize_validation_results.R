# This is a Rscript comparing perturbPlan results with QC downsampling results only
# Modified version with light/heavy tail only and simplified formatting
library(dplyr) 
library(tidyr) 
library(ggplot2)
library(patchwork)
library(stringr)
library(cowplot)
library(scales)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
experiment <- args[1]

# specify significance level
if(experiment == "Gasperini"){
  alpha <- 0.01
  window_def <- "gene_analysis"
}else{
  alpha <- 0.1
  window_def <- "enhancer_analysis"
}

# source the necessary scripts
source("~/.Rprofile")
source("reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R")

# load the simulation results
intermediate_files_folder <- "reproduction-code-prep/realdata-validation-pipeline/intermediate-files"

########################## 0. Prepare for summarizing ##########################
# input the parameter dataframe
parameter_df <- readRDS(sprintf("%s/%s/real_data_parameter_grid.rds", intermediate_files_folder, experiment))
result_dir <- paste0(.get_config_path("LOCAL_PERTURBPLAN_DATA_DIR"), "/realdata-validation-working")
simulation_results <- readRDS(sprintf("%s/%s/real_data_results.rds", result_dir, experiment))

# load the analytic results
analytic_results <- readRDS(sprintf("reproduction-code-prep/realdata-validation-pipeline/PerturbPlan/results/%s_%s_perturbplan_overall_power.rds",
                                    experiment, window_def))

######################## 1. Perform multiplicity correction ####################
# perform multiplicity correction with memory optimization
simulation_power <- NULL
fc_tail_list <- c("light_tail", "heavy_tail")  # Only light and heavy tail

# Pre-process simulation results once to avoid repeated unnesting
cat("Pre-processing simulation results...\n")
processed_sim_results <- simulation_results$results |>
  tidyr::unnest_longer(output) |>
  tidyr::unnest(output) |>
  dplyr::mutate(pair_id = names(output)) |>
  dplyr::rename(p_value = output, QC_type = output_id)

# Clear original large object
rm(simulation_results)
gc()

# Helper function to compute simulation power for QC type only
compute_simulation_power_qc <- function() {
  simulation_power_temp <- NULL
  
  for (fc_tail in fc_tail_list) {
    cat("Processing", fc_tail, "tail for downsampling with QC...\n")
    
    # load the discovery pairs
    discovery_set <- readRDS(sprintf("%s/%s/discovery_pairs_%s_%s.rds",
                                     intermediate_files_folder, experiment, fc_tail, window_def))

    # compute simulation power for this tail type (QC only)
    tail_power <- processed_sim_results |>
      dplyr::filter(pair_id %in% discovery_set$pair_id) |>
      dplyr::group_by(run_id, grid_id, QC_type) |>
      dplyr::mutate(
        rejection = tidyr::replace_na((p.adjust(p_value, method = "BH") <= alpha), FALSE),
        QC = is.na(p_value)
      ) |>
      dplyr::ungroup() |>
      dplyr::group_by(pair_id, grid_id, QC_type) |>
      dplyr::summarise(
        sceptre_power = mean(rejection, na.rm = TRUE),
        sceptre_QC = mean(QC),
        .groups = 'drop'
      ) |>
      dplyr::left_join(parameter_df |> dplyr::select(R, N, grid_id), by = "grid_id") |>
      dplyr::left_join(discovery_set, by = "pair_id") |>
      dplyr::rename(num_total_cells = N, reads_per_cell = R) |> 
      dplyr::filter(pair_type == "positive") |>
      dplyr::filter(QC_type == "p_value_w_QC") |>  # Only QC version
      dplyr::group_by(num_total_cells, reads_per_cell) |> 
      dplyr::summarise(simulation_power = mean(sceptre_power), .groups = 'drop') |>
      dplyr::mutate(tail_type = fc_tail, method = "downsampling_with_qc")
    
    # Bind results
    simulation_power_temp <- dplyr::bind_rows(simulation_power_temp, tail_power)
    rm(tail_power, discovery_set)
  }
  return(simulation_power_temp)
}

# Compute simulation power for QC only
simulation_power <- compute_simulation_power_qc()
gc()

########################### 2. Merge dataframes with analytic approaches ############################
# Combine both lower bound and median estimate analytic results (filter for light/heavy only)
combined_analytic_results <- analytic_results |>
  dplyr::filter(tail_type %in% fc_tail_list) |>  # Only light and heavy tail
  dplyr::select(reads_per_cell, num_total_cells, num_trt_cells, tail_type, power_type, power) |>
  tidyr::pivot_wider(names_from = power_type, values_from = power) |>
  tidyr::pivot_longer(cols = c("lower_bound", "median_estimate"), 
                      names_to = "power_type", values_to = "analytic_power")

# Create merged results with both power types and analytic approaches
merged_results <- combined_analytic_results |>
  dplyr::inner_join(simulation_power, by = c("tail_type", "num_total_cells", "reads_per_cell")) |>
  dplyr::select(num_total_cells, num_trt_cells, reads_per_cell, tail_type, 
                power_type, analytic_power, simulation_power, method) |>
  dplyr::mutate(
    reads_per_cell = factor(
      as.character(reads_per_cell),
      levels = as.character(sort(unique(reads_per_cell)))
    ),
    num_total_cells = factor(
      as.character(num_total_cells),
      levels = as.character(sort(unique(num_total_cells)))
    ),
    num_trt_cells = factor(
      as.character(num_trt_cells),
      levels = as.character(sort(unique(num_trt_cells)))
    ),
    tail_type = ifelse(tail_type == "light_tail", "Lower effect variability", "Higher effect variability"),
    tail_type = factor(tail_type, levels = c("Lower effect variability", "Higher effect variability"))
  ) |> 
  # Rename and combine power columns
  dplyr::rename(downsampling_power = simulation_power) |>
  tidyr::pivot_longer(cols = c("analytic_power", "downsampling_power"),
                      names_to = "result_type", values_to = "power") |>
  # Create unified method column
  dplyr::mutate(
    power_method = case_when(
      result_type == "analytic_power" & power_type == "lower_bound" ~ "PerturbPlan (min_FC)",
      result_type == "analytic_power" & power_type == "median_estimate" ~ "PerturbPlan (median_fc)",
      result_type == "downsampling_power" ~ "Downsampling",  # Simplified name
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(power_method)) |>
  dplyr::select(-result_type, -power_type, -method)

# compute the error metric
range(merged_results |> dplyr::filter(power_method != "PerturbPlan (min_FC)") |> 
        dplyr::distinct() |> 
        tidyr::pivot_wider(values_from = "power", names_from = "power_method", 
                           id_cols = c("reads_per_cell", "num_total_cells", "tail_type")) |> 
        dplyr::mutate(abs_error = abs(Downsampling - `PerturbPlan (median_fc)`)) |> 
        dplyr::select(abs_error) |> dplyr::pull())

# Absolute error
# Gasperini: 0.001050032 0.119722895
# Ray: 0.0008653093 0.1189883860


############################# 3. Visualize the results with facet_grid #########################
# Get default values from parameter_df where arm = FALSE
default_total_cells <- parameter_df |>
  dplyr::filter(arm_N == FALSE) |>
  dplyr::pull(N) |>
  unique() |>
  as.character()

default_reads_per_cell <- parameter_df |>
  dplyr::filter(arm_R == FALSE) |>
  dplyr::pull(R) |>
  unique() |>
  as.character()

# Create combined plot with two rows: one for total cells, one for reads per cell
# First row: facet by cells per target (fixing reads per cell)
by_cells_per_target <- merged_results |>
  dplyr::filter(reads_per_cell == default_reads_per_cell) |>
  ggplot(aes(x = num_trt_cells, y = power, color = power_method, group = power_method)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ tail_type, scales = "free_x", 
             labeller = labeller(tail_type = function(x) str_replace_all(str_to_title(x), "_", " "))) +  # Simplified strip text
  labs(
    x = "Cells per target",
    y = "Power",
    color = "Method"
  ) +
  scale_color_manual(
    values = c("PerturbPlan (min_FC)" = "#fee090", 
               "PerturbPlan (median_fc)" = "#d73027",
               "Downsampling" = "#74add1"),  # Simplified legend
    labels = c("PerturbPlan (min_FC)" = "PerturbPlan (min FC)", 
               "PerturbPlan (median_fc)" = "PerturbPlan (median FC)",
               "Downsampling" = "Downsampling")  # Simplified legend
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  )

# Second row: facet by reads per cell (fixing number of total cells)
by_reads_per_cell <- merged_results |>
  dplyr::filter(num_total_cells == default_total_cells) |>
  ggplot(aes(x = reads_per_cell, y = power, color = power_method, group = power_method)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ tail_type, scales = "free_x", 
             labeller = labeller(tail_type = function(x) str_replace_all(str_to_title(x), "_", " "))) +  # Simplified strip text
  labs(
    x = "Reads per cell",
    y = "Power",
    color = "Method"
  ) +
  scale_color_manual(
    values = c("PerturbPlan (min_FC)" = "#fee090", 
               "PerturbPlan (median_fc)" = "#d73027",
               "Downsampling" = "#74add1"),  # Simplified legend
    labels = c("PerturbPlan (min_FC)" = "PerturbPlan (min FC)", 
               "PerturbPlan (median_fc)" = "PerturbPlan (median FC)",
               "Downsampling" = "Downsampling")  # Simplified legend
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  )

# Third row: effect size distribution for each discovery set (light and heavy only)
# Load the discovery data for each tail type to get the actual effect sizes used
light_tail_data <- readRDS(sprintf("reproduction-code-prep/realdata-validation-pipeline/intermediate-files/%s/fc_distr_light_tail_%s.rds", experiment, window_def))
heavy_tail_data <- readRDS(sprintf("reproduction-code-prep/realdata-validation-pipeline/intermediate-files/%s/fc_distr_heavy_tail_%s.rds", experiment, window_def))

# Get the specific plugin effect sizes used for analytic power calculation
light_max_effect <- max(light_tail_data$fold_change, na.rm = TRUE)
light_median_effect <- median(light_tail_data$fold_change, na.rm = TRUE)
heavy_max_effect <- max(heavy_tail_data$fold_change, na.rm = TRUE)
heavy_median_effect <- median(heavy_tail_data$fold_change, na.rm = TRUE)

# Calculate common x-axis range across light and heavy discovery sets
all_fold_changes <- c(light_tail_data$fold_change, heavy_tail_data$fold_change)
x_range <- range(all_fold_changes, na.rm = TRUE)
x_limits <- c(x_range[1] - 0.01, x_range[2] + 0.01)  # Add small buffer

# Create fold change data with tail type for faceting
fold_change_data <- rbind(
  data.frame(fold_change = light_tail_data$fold_change, 
             tail_type = "light_tail",
             max_effect = light_max_effect,
             median_effect = light_median_effect),
  data.frame(fold_change = heavy_tail_data$fold_change, 
             tail_type = "heavy_tail",
             max_effect = heavy_max_effect,
             median_effect = heavy_median_effect)
) |>
  dplyr::mutate( 
    tail_type = ifelse(tail_type == "light_tail", "Lower effect variability", "Higher effect variability"),
    tail_type = factor(tail_type, levels = c("Lower effect variability", "Higher effect variability")))

# Create single faceted plot using same approach as power plots
effect_size_plot <- fold_change_data |>
  dplyr::mutate(fold_change = (1 - fold_change) * 100) |>
  ggplot(aes(x = fold_change)) +
  xlim(sort((1 - x_limits) * 100)) +
  geom_histogram(bins = 30, alpha = 0.7, fill = "lightgrey", color = "black") +
  # Add lines for both max and median effect sizes
  geom_vline(aes(xintercept = (1 - max_effect) * 100), 
             color = "#fee090", linetype = "solid", linewidth = 1.2) +
  geom_vline(aes(xintercept = (1 - median_effect) * 100), 
             color = "#d73027", linetype = "solid", linewidth = 1.2) +
  facet_wrap(~ tail_type, scales = "free_x",
             labeller = labeller(tail_type = function(x) str_replace_all(str_to_title(x), "_", " "))) +  # Same as power plots
  labs(
    x = "Fold Change (percent)",
    y = "Count",
    title = ifelse(experiment == "Gasperini", "Perturb-seq validation", "TAP-seq validation")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14)
  )

# Extract legend from the power plot for separate use
power_plot_with_legend <- by_reads_per_cell  # This has the legend
legend_plot <- get_plot_component(power_plot_with_legend + theme(legend.title = element_blank()), 
                                  "guide-box", return_all = TRUE)[[3]]

# Combine all three plots vertically with effect size plot at the top (no title)
validation_plot <- effect_size_plot / by_cells_per_target / by_reads_per_cell

# Save the ggplot object as RDS to results/ folder
results_path <- "reproduction-code-prep/realdata-validation-pipeline/PerturbPlan/results"
saveRDS(validation_plot, sprintf("%s/%s_%s_plot.rds", results_path, experiment, window_def))

# Save the legend as separate RDS file
saveRDS(legend_plot, sprintf("%s/%s_%s_legend_plot.rds", results_path, experiment, window_def))
