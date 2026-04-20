# This is a script investigating the optimal design fixing cost or power or both
library(metR)
library(scales)
library(ggplot2)
library(patchwork)
library(cowplot)

############################# Plot for varying TPM #############################
# constract label dataframe
label_df <- cost_grid_opt_tpm |>                           # one point per (assay, cost)
  group_by(cost_of_interest) |>             
  slice_max(cells_per_target, n = 1, with_ties = FALSE) %>% 
  ungroup()

# plot the Perturb-seq
opt_tpm_plot <- cost_power_varying_tpm_df |> 
  ggplot(aes(x = cells_per_target, y = sequenced_reads_per_cell, color = TPM_threshold)) +
  geom_smooth(se = FALSE) +
  geom_smooth(
    data        = cost_grid_opt_tpm,
    mapping     = aes(x = cells_per_target,          # ← add x
                      y = sequenced_reads_per_cell,            # ← add y
                      linetype = as.factor(cost_of_interest)),
    se = FALSE,
    color = "black",
    inherit.aes = FALSE
  ) +
  ## labels on the curves
  geom_text(data  = label_df,
            aes(x = cells_per_target,     
                y = sequenced_reads_per_cell,            
                label = dollar(cost_of_interest),
                group  = cost_of_interest),
            hjust = 0.4,           # nudge a little to the right of the line end
            vjust = 0.8,
            size  = 3,
            inherit.aes = FALSE) +
  ## linetype legend title
  scale_linetype_discrete(name = "Total cost ($)") +
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = "Cells per target", y = "Reads per cell", color = "TPM analysis threshold",
       title = sprintf("Optimal TPM threshold given power at %.2f", power_target)) +
  geom_vline(xintercept = cells_per_target, color = "#00BFC4", linetype = "dashed") +
  geom_hline(yintercept = sequenced_reads_per_cell, color = "#00BFC4", linetype = "dashed") +
  theme_bw() +
  ## remove the linetype legend, keep the colour one
  scale_linetype_discrete(guide = "none") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position   = "bottom",
    plot.title = element_text(hjust = 0.5, size = 20)
  )

# obtain optimal cost df
optimal_cost_df <- cost_power_varying_tpm_df |>
  dplyr::filter(total_cost == minimum_cost) |>
  dplyr::mutate(TPM_threshold = as.numeric(as.character(TPM_threshold)))

# plot the reads over effect size
read_tpm_plot <- optimal_cost_df |>
  ggplot(aes(x = TPM_threshold, y = sequenced_reads_per_cell)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "TPM analysis threshold", y = "Reads per cell") +
  geom_hline(yintercept = sequenced_reads_per_cell, color = "#00BFC4", linetype = "dashed") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position   = "bottom"
  )
cell_tpm_plot <- optimal_cost_df |>
  ggplot(aes(x = TPM_threshold, y = cells_per_target)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "TPM analysis threshold", y = "Cells per target") +
  geom_hline(yintercept = cells_per_target, color = "#00BFC4", linetype = "dashed") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position   = "bottom"
  )
cost_tpm_plot <- optimal_cost_df |>
  ggplot(aes(x = TPM_threshold, y = minimum_cost)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "TPM analysis threshold", y = "Optimal cost") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position   = "bottom"
  )

################################# Plot for varying FC ##########################
# construct label dataframe
label_df <- cost_grid_opt_fc |>                           # one point per (assay, cost)
  group_by(cost_of_interest) |>             # ──► adjust if you have both assays
  slice_max(cells_per_target, n = 1, with_ties = FALSE) %>% 
  ungroup()

# plot the Perturb-seq
opt_fc_plot <- cost_power_varying_fc_df |> 
  ggplot(aes(x = cells_per_target, y = sequenced_reads_per_cell, color = minimum_fold_change)) +
  geom_smooth(se = FALSE) +
  geom_smooth(
    data        = cost_grid_opt_fc,
    mapping     = aes(x = cells_per_target,          # ← add x
                      y = sequenced_reads_per_cell,            # ← add y
                      linetype = as.factor(cost_of_interest)),
    se = FALSE,
    color = "black",
    inherit.aes = FALSE
  ) +
  ## labels on the curves
  geom_text(data  = label_df,
            aes(x = cells_per_target,     
                y = sequenced_reads_per_cell,            
                label = dollar(cost_of_interest),
                group  = cost_of_interest),
            hjust = 0.5,           # nudge a little to the right of the line end
            vjust = 0.8,
            size  = 3,
            inherit.aes = FALSE) +
  ## linetype legend title
  scale_linetype_discrete(name = "Total cost ($)") +
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = "Cells per target", y = "Reads per cell", color = "Minimum effect size (percent)", 
       title = sprintf("Optimal effect size given power at %.2f", power_target)) +
  geom_vline(xintercept = cells_per_target, color = "#00BFC4", linetype = "dashed") +
  geom_hline(yintercept = sequenced_reads_per_cell, color = "#00BFC4", linetype = "dashed") +
  theme_bw() +
  ## remove the linetype legend, keep the colour one
  scale_linetype_discrete(guide = "none") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position   = "bottom",
    plot.title = element_text(hjust = 0.5, size = 20)
  )

# obtain optimal cost dataframe
optimal_cost_df <- cost_power_varying_fc_df |>
  dplyr::filter(total_cost == minimum_cost) |>
  dplyr::mutate(minimum_fold_change = as.numeric(as.character(minimum_fold_change)))

# plot the reads over effect size
read_fc_plot <- optimal_cost_df |>
  ggplot(aes(x = minimum_fold_change, y = sequenced_reads_per_cell)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "Minimum effect size (percent)", y = "Reads per cell") +
  geom_hline(yintercept = sequenced_reads_per_cell, color = "#00BFC4", linetype = "dashed") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
cell_fc_plot <- optimal_cost_df |>
  ggplot(aes(x = minimum_fold_change, y = cells_per_target)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "Minimum effect size (perccent)", y = "Cells per target") +
  geom_hline(yintercept = cells_per_target, color = "#00BFC4", linetype = "dashed") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
cost_fc_plot <- optimal_cost_df |>
  ggplot(aes(x = minimum_fold_change, y = minimum_cost)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "Minimum effect size (percent)", y = "Optimal cost") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

# combined two plots
combined_opt_tpm_fc_plot <- ((opt_tpm_plot + opt_fc_plot +  plot_layout(guides = "collect")) &
                               theme(legend.position = "bottom"))

############################### Optimal cost over TPM/FC #######################
# combined four plots varying TPM/minimum FC
combined_opt_cost_cell_read_plot <- ((cell_tpm_plot + cell_fc_plot) / (read_tpm_plot + read_fc_plot) / (cost_tpm_plot + cost_fc_plot) )  +
  plot_annotation(
    title = "Cells, reads and cost over TPM threshold or effect size given power at optimal cost",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 20))
  )
