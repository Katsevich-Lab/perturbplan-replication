# This is a script comparing the simulated power and anlytic power
library(ggplot2)
library(dplyr)
library(perturbplan)
library(patchwork)

# load the baseline information
baseline_expression <- readRDS("reproduction-code-prep/simulation-benchmarking/baseline_expression.rds")
simulated_power <- readRDS("reproduction-code-prep/simulation-benchmarking/averaged_simulation_results.rds")

# prepare for analytic power computation 
analytic_power <- NULL
timing_data <- NULL

for (param in 1:nrow(simulated_power)) {
  
  # extract the parameter grid
  library_size <- sum(baseline_expression$mean)
  gRNAs_per_target <- 15
  multiple_testing_alpha <- 0.1
  multiple_testing_method <- "BH"
  control_group <- "complement"
  side <- "left"
  fold_change_mean <- 1 - simulated_power$simulation_effect_size[param]
  fold_change_sd <- 0.13
  prop_non_null <- 0.05
  
  # extract the fc_expression_df
  custom_expression <- list(
    baseline_expression = list(baseline_expression = baseline_expression |> 
                                 dplyr::mutate(relative_expression = mean / library_size,
                                               expression_size = 1 / dispersion) |> 
                                 dplyr::rename(response_id = gene_id) |>
                                 dplyr::select(response_id, relative_expression, expression_size))
  )
  fc_expression_df <- extract_fc_expression_info(minimum_fold_change = fold_change_mean,
                                                 gRNA_variability = fold_change_sd,
                                                 B = 1000, TPM_threshold = 0,
                                                 custom_pilot_data = custom_expression,
                                                 gRNAs_per_target = gRNAs_per_target)
  
  # compute the analytical overall power with timing
  num_trt_cells <- simulated_power$num_treatment_cells[param]
  num_cntrl_cells <- simulated_power$num_control_cells[param]
  
  # Record timing for PerturbPlan computation
  start_time <- Sys.time()
  power_result <- compute_power_plan_overall(
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    library_size = library_size,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    side = side,
    fc_expression_df = fc_expression_df$fc_expression_df,
    prop_non_null = prop_non_null,
    return_full_results = TRUE
  )
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Store timing data
  timing_record <- data.frame(
    effect_size = 1 - fold_change_mean,
    num_treatment_cells = num_trt_cells,
    start_time = start_time,
    end_time = end_time,
    duration_seconds = duration,
    duration_minutes = duration / 60,
    duration_hours = duration / 3600
  )
  timing_data <- rbind(timing_data, timing_record)
  
  analytic_power <- data.frame(power_result) |> 
    dplyr::rename(num_treatment_cells = num_trt_cells) |>
    dplyr::mutate(effect_size = 1 - fold_change_mean) |> 
    dplyr::bind_rows(analytic_power)
}

# concatenate the results - round effect_size to 2 digits for alignment
merged_results <- simulated_power |> 
  dplyr::select(num_treatment_cells, simulation_effect_size, averaged_power) |>
  dplyr::mutate(effect_size = round(simulation_effect_size, 2)) |>
  dplyr::left_join(analytic_power |> 
                     dplyr::select(overall_power, num_treatment_cells, effect_size) |>
                     dplyr::mutate(effect_size = round(effect_size, 2)) |>
                     dplyr::distinct(), 
                   by = c("num_treatment_cells", "effect_size")) |> 
  dplyr::rename(analytic_power = overall_power,
                simulated_power = averaged_power) |>
  dplyr::select(num_treatment_cells, effect_size, simulated_power, analytic_power)

# reshape data using pivot_longer
plot_data <- merged_results |>
  tidyr::pivot_longer(cols = c("simulated_power", "analytic_power"),
                      names_to = "method", 
                      values_to = "power") |>
  dplyr::mutate(fold_change = 1 - effect_size)  # Calculate fold change

# Debug: Check what effect_size values are near 0.15
cat("Effect sizes near 0.15:\n")
debug_data <- plot_data |> dplyr::filter(abs(effect_size - 0.15) < 0.01)
print(unique(debug_data$effect_size))

# First plot: Power vs Cells per target (fixing effect_size = 0.15)
# Use the exact effect_size closest to 0.15
target_effect_size <- plot_data |> 
  dplyr::filter(abs(effect_size - 0.15) < 0.01) |>
  dplyr::pull(effect_size) |>
  unique() |>
  {\(x) x[which.min(abs(x - 0.15))]}()

cat("Using effect_size =", target_effect_size, "\n")

plot1_data <- plot_data |>
  dplyr::filter(effect_size == target_effect_size) |>  # Use exact match
  dplyr::group_by(num_treatment_cells, method) |>
  dplyr::summarise(power = mean(power, na.rm = TRUE), .groups = 'drop')

power_vs_cells <- ggplot(plot1_data, aes(x = num_treatment_cells, y = power, color = method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line(aes(group = method), linewidth = 1) +
  labs(
    x = "Cells per target",
    y = "Power",
    color = "Method"
  ) +
  scale_color_manual(
    values = c("analytic_power" = "#d73027", "simulated_power" = "#4575b4"),
    labels = c("analytic_power" = "PerturbPlan", "simulated_power" = "Simulation")
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.75, 0.3),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

# Second plot: Power vs Fold change (fixing num_treatment_cells = 1000)
plot2_data <- plot_data |>
  dplyr::mutate(fold_change = round((1 - fold_change) * 100)) |>
  dplyr::filter(abs(num_treatment_cells - 1000) < 50) |>  # Filter for cells ≈ 1000
  dplyr::distinct(effect_size, method, .keep_all = TRUE) |>  # Remove exact duplicates
  dplyr::group_by(effect_size, method) |>
  dplyr::summarise(power = mean(power, na.rm = TRUE), 
                   fold_change = mean(fold_change, na.rm = TRUE),
                   .groups = 'drop')

# compute absolute error
max(
  max(plot1_data |> 
        tidyr::pivot_wider(values_from = "power", id_cols = "num_treatment_cells", names_from = "method") |> 
        dplyr::mutate(abs_error = abs(analytic_power - simulated_power)) |> 
        dplyr::select(abs_error) |> dplyr::pull()),
  max(plot2_data |> tidyr::pivot_wider(values_from = "power", id_cols = "fold_change", names_from = "method") |> 
        dplyr::mutate(abs_error = abs(analytic_power - simulated_power)) |> 
        dplyr::select(abs_error) |> dplyr::pull())
)

power_vs_fold_change <- ggplot(plot2_data, aes(x = fold_change, y = power, color = method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line(aes(group = method), linewidth = 1) +
  labs(
    x = "Fold change (percent)",
    y = "Power",
    color = "Method", 
    title = "Simulation benchmarking"
  ) +
  scale_color_manual(
    values = c("analytic_power" = "#d73027", "simulated_power" = "#4575b4"),
    labels = c("analytic_power" = "PerturbPlan", "simulated_power" = "Simulation")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.text = element_text(size = 12)
  )

# Third plot: Timing comparison (averaged)
# Load simulation timing data
simulation_timing <- readRDS("code/benchmarking/simulation_timing_data.rds")

# Calculate average timing for each method (scaled by 1000)
avg_timing <- data.frame(
  method = c("PerturbPlan", "Simulation"),
  avg_duration_scaled = c(
    mean(timing_data$duration_seconds, na.rm = TRUE) * 1000,
    mean(simulation_timing$total_seconds, na.rm = TRUE) * 1000
  )
) |>
  dplyr::mutate(method = factor(method, levels = c("PerturbPlan", "Simulation")))

# Create timing comparison plot with averaged values
timing_plot <- ggplot(avg_timing, aes(x = method, y = avg_duration_scaled, fill = method)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(
    x = "Method",
    y = "Runtime (seconds)",
    fill = "Method"
  ) +
  scale_fill_manual(
    values = c("PerturbPlan" = "#d73027", "Simulation" = "#4575b4")
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Will share legend with other plots
  )

# First combine the two power plots with shared legend
power_plots <- (power_vs_cells | power_vs_fold_change) 

# Then add the timing plot without legend
combined_plot <- power_plots | timing_plot

# Save the plotting object
saveRDS(combined_plot, "reproduction-code-prep/simulation-benchmarking/benchmarking_plot.rds")
