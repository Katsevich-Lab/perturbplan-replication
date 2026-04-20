# Reproduce PerturbPlan "Cost Minimization" Plot
# This script uses the actual perturbplan package to generate real cost-power tradeoff curves
# Scenario: Minimize total cost while varying cells per target and reads per cell,
# with power constraint = 0.8

library(perturbplan)
library(ggplot2)
library(scales)
library(stringr)
library(grid)

# Extract expression info (using K562 data)
expression_info <- perturbplan:::extract_expression_info(
  biological_system = "K562",
  B = 1000,              # Sample 1000 genes
  gene_list = NULL,      # Random sampling
  TPM_threshold = 1,     # TPM filtering
  custom_pilot_data = NULL
)

# Set up parameters for cost minimization
# Both cells_per_target and reads_per_cell will vary to minimize cost
params <- list(
  minimizing_variable = "cost",  # Minimize cost
  fixed_variable = list(
    TPM_threshold = 1,               # Fixed at default value
    minimum_fold_change = 0.8        # Fixed at default value
    # cells_per_target and reads_per_cell are NOT in fixed_variable - they vary
  ),
  
  # Experimental parameters (defaults)
  MOI = 10,
  num_targets = 1000,
  non_targeting_gRNAs = 50,
  gRNAs_per_target = 4,
  
  # Effect size parameters (defaults)
  gRNA_variability = 0,
  prop_non_null = 0.25,
  
  # Analysis parameters (defaults)
  control_group = "complement",
  side = "both",
  multiple_testing_alpha = 0.1,
  
  # Power constraint
  power_target = 0.8,
  
  # Cost parameters
  cost_constraint = NULL,  # No cost constraint - we're minimizing cost
  cost_per_captured_cell = 0.086,
  cost_per_million_reads = 1,
  
  # Grid parameters
  grid_size = 100,  # Larger grid for 2D optimization
  
  # Mapping efficiency
  mapping_efficiency = 0.5,
  
  # Pilot data
  baseline_expression_stats = expression_info$expression_df,
  library_parameters = expression_info$pilot_data$library_parameters
)

# Step 1: Call cost_power_computation to get power-cost grid
cat("Running perturbplan cost-power computation...\n")
cost_power_grid <- do.call(perturbplan:::cost_power_computation, params)

# Step 2: Call find_optimal_cost_design to generate equi-power and equi-cost curves
cat("Finding optimal cost design...\n")
optimal_results <- perturbplan:::find_optimal_cost_design(
  cost_power_df = cost_power_grid,
  minimizing_variable = "cost",
  power_target = 0.8,
  power_precision = 0.002,
  MOI = params$MOI,
  num_targets = params$num_targets,
  non_targeting_gRNAs = params$non_targeting_gRNAs,
  gRNAs_per_target = params$gRNAs_per_target,
  cost_per_captured_cell = params$cost_per_captured_cell,
  cost_per_million_reads = params$cost_per_million_reads,
  cost_grid_size = 100
)

# Extract data for plotting
power_data <- optimal_results$optimal_cost_power_df |> 
  dplyr::filter(cells_per_target <= 2500) # Equi-power curves
cost_data <- optimal_results$optimal_cost_grid        # Equi-cost curves

# Find optimal point (minimum cost at target power)
target_power <- 0.8
power_tolerance <- 0.01
target_rows <- power_data[abs(power_data$overall_power - target_power) <= power_tolerance, ]
optimal_idx <- which.min(target_rows$total_cost)
optimal_point <- target_rows[optimal_idx, ]

cat("\n=== Optimal Solution ===\n")
cat("Cells per target:", comma(round(optimal_point$cells_per_target)), "\n")
cat("Reads per cell:", comma(round(optimal_point$sequenced_reads_per_cell)), "\n")
cat("Total cost: $", comma(round(optimal_point$total_cost, 2)), "\n")
cat("Achieved power:", percent(optimal_point$overall_power, accuracy = 0.1), "\n")

# Format optimal cost for legend (round to nearest thousand)
optimal_cost_k <- round(optimal_point$total_cost / 1000)
cost_legend_label <- paste0("Optimal cost: $", optimal_cost_k, "k")

# Define colors
equi_power_color <- "#2E86AB"  # Blue for equi-power
equi_cost_color <- "#A23B72"   # Purple-red for equi-cost

# (Optional) sanitize label in case it has stray spaces
cost_legend_label <- str_trim(cost_legend_label)

p <- ggplot() +
  # Equi-power (solid) as a color legend entry
  geom_smooth(
    data = power_data,
    aes(x = cells_per_target, y = sequenced_reads_per_cell, color = "Power \u2265 0.8"),
    se = FALSE, size = 0.8
  ) +
  geom_point(
    data = power_data,
    aes(x = cells_per_target, y = sequenced_reads_per_cell, color = "Power \u2265 0.8"),
    size = 1
  ) +
  
  # Equi-cost (dashed) as a color legend entry
  geom_smooth(
    data = cost_data,
    aes(x = cells_per_target, y = sequenced_reads_per_cell, color = cost_legend_label),
    se = FALSE, size = 0.8
  ) +
  geom_point(
    data = cost_data,
    aes(x = cells_per_target, y = sequenced_reads_per_cell, color = cost_legend_label),
    size = 1
  ) +
  
  # Optimal point as a shape legend entry
  geom_point(
    data = optimal_point,
    aes(x = cells_per_target, y = sequenced_reads_per_cell, shape = "Optimal design"),
    color = "darkgoldenrod", size = 4
  ) +
  
  # Color legend for power and cost
  scale_color_manual(
    name = "",
    values = setNames(
      c(equi_power_color, equi_cost_color),
      c("Power \u2265 0.8", cost_legend_label)
    ),
    breaks = c("Power \u2265 0.8", cost_legend_label)
  ) +
  
  # Shape legend for optimal point
  scale_shape_manual(
    name = "",
    values = c("Optimal design" = 18),
    guide = guide_legend(order = 2)
  ) +
  
  # Merge the two guides area + tidy the color legend keys
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 1.2, shape = NA)),
    shape = guide_legend(order = 2)
  ) +
  
  # Axes
  scale_y_log10(labels = c("30k", "100k", "300k"), breaks = c(30000, 100000, 300000)) +

  labs(x = "Cells per target", y = "Reads per cell",
       title = "Minimize total cost") +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.box = "vertical",
    legend.margin = margin(t = 0, b = 0),
    #    legend.box.margin = margin(t = -3, b = -3),
    legend.title = element_blank(),
    legend.spacing.y  = unit(0, "pt"),
    #    legend.key.spacing.y = unit(-2, "pt"),
    legend.text  = element_text(margin = margin(b = 0)),
    plot.margin = unit(c(5, 15, 5, 5), "pt")  # bump up right margin
  )

p

# Save the plot
ggsave("reproduction-code-prep/final_plots/figures/cost_minimization_plot.pdf", p, width = 3.4, height = 2.5)
