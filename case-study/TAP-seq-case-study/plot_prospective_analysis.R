# This is a script investigating the optimal design fixing cost or power or both
library(metR)
library(scales)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrepel)

############### Cells, reads and cost over primer efficiency ###################
# load parameter grid
parameter_grid <- readRDS(path_to_param_df)

# obtain the optimal cost dataframe
optimal_cost_eff_df <- optimal_power_df_eff |>
  dplyr::group_by(assay, primer_efficiency_bound, cost_ratio) |>
  dplyr::filter(total_cost == min(total_cost)) |>
  dplyr::ungroup() |>
  dplyr::mutate(minimum_fold_change = as.numeric(as.character(minimum_fold_change))) |> 
  dplyr::select(cost_ratio, total_cost, assay, primer_efficiency_bound, cells_per_target, sequenced_reads_per_cell) |>
  dplyr::mutate(primer_filtering = primer_efficiency_bound) |>
  dplyr::left_join(parameter_grid |> dplyr::rename(cost_ratio = ratio), 
                   by = c("cost_ratio", "primer_efficiency_bound"))

# plot the reads, cells and cost over primer filtering
read_eff_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_primer_efficiency_bound) |>
  ggplot(aes(x = primer_filtering, y = sequenced_reads_per_cell, color = assay)) +
  geom_point() +
  scale_y_log10() +
  geom_line(linewidth = 1) +
  labs(x = "Primer efficiency", y = "Reads per cell",
       title = "Optimal reads versus primer efficiency") +
  geom_hline(yintercept = reads_per_cell, color = "#00BFC4", linetype = "dashed") +
  annotate(
    "text",
    x      = 0.1,   # 80 % of the way across
    y      = reads_per_cell,
    label  = sprintf("Reads per cell = %d", round(reads_per_cell)),
    hjust  = 0.6, vjust = 1.2,    # tweak so it sits just above the line
    size   = 4, colour = "#00BFC4"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16))
cell_eff_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_primer_efficiency_bound) |>
  ggplot(aes(x = primer_filtering, y = cells_per_target, color = assay)) +
  geom_point() +
  geom_line(linewidth = 1) +
  labs(x = "Primer efficiency", y = "Cells per target",
       title = "Optimal cells versus primer efficiency") +
  geom_hline(yintercept = cells_per_target, color = "#00BFC4", linetype = "dashed") +
  annotate(
    "text",
    x      = 0.1,
    y      = cells_per_target,
    label  = sprintf("Cells per target = %d", round(cells_per_target)),
    hjust  = 0.6, vjust = -0.3,    # tweak so it sits just above the line
    size   = 4, colour = "#00BFC4"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16))
cost_eff_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_primer_efficiency_bound) |>
  ggplot(aes(x = primer_filtering, y = total_cost, color = assay)) +
  geom_point() +
  geom_line(linewidth = 1) +
  labs(x = "Primer efficiency", y = "Optimal cost ($)",
       title = "Optimal cost versus primer efficiency") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  )


################ Cells, reads and cost over cost ratio plot ####################
# plot the reads, cells and cost over cost ratio
read_ratio_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_ratio) |>
  ggplot(aes(x = cost_ratio, y = sequenced_reads_per_cell, color = assay)) +
  geom_point() +
  scale_x_log10() +
  geom_line(linewidth = 1) +
  labs(x = "Library cost / Sequencing cost", y = "Reads per cell",
       title = "Optimal reads versus cost ratio") +
  geom_hline(yintercept = reads_per_cell, color = "#00BFC4", linetype = "dashed") +
  annotate(
    "text",
    x      = 0.5,   
    y      = reads_per_cell,
    label  = sprintf("Reads per cell = %d", round(reads_per_cell)),
    hjust  = 0.6, vjust = 1.2,    # tweak so it sits just above the line
    size   = 4, colour = "#00BFC4"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16))
cell_ratio_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_ratio) |>
  ggplot(aes(x = cost_ratio, y = cells_per_target, color = assay)) +
  geom_point() +
  scale_x_log10() +
  geom_line(linewidth = 1) +
  labs(x = "Library cost / Sequencing cost", y = "Cells per target",
       title = "Optimal cells versus cost ratio") +
  geom_hline(yintercept = cells_per_target, color = "#00BFC4", linetype = "dashed") +
  annotate(
    "text",
    x      = 0.6,
    y      = cells_per_target,
    label  = sprintf("Cells per target = %d", round(cells_per_target)),
    hjust  = 0.6, vjust = -0.3,    # tweak so it sits just above the line
    size   = 4, colour = "#00BFC4"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16))
cost_ratio_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_ratio) |>
  ggplot(aes(x = cost_ratio, y = total_cost, color = assay)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_line(linewidth = 1) +
  labs(x = "Library cost / Sequencing cost", y = "Optimal cost ($)",
       title = "Optimal cost versus cost ratio") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

###################### Optimal cost over efficiency/cost ratio #################
## 1 ── Re-make the three panels but turn their legends off
panels <- (((cell_eff_plot + cell_ratio_plot) / (read_eff_plot + read_ratio_plot) / (cost_eff_plot + cost_ratio_plot) ) & 
             theme(legend.position = "none")) +
  plot_annotation(
    title = sprintf("Cells, reads and cost over primer efficiency or cost ratio given power %.2f at optimal cost", power_target),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

## 2 ── Stack the panels over the shared legend
combined_opt_cost_cell_read_plot_over_eff_and_cost_ratio <- panels / shared_legend + plot_layout(heights = c(0.3, 0.3, 0.3, 0.05)) 

# save plotting object
saveRDS(panels, paste0(dir_to_results, "combined_cost_read_cell_eff_cost_ratio_plotting.rds"))

########################## Plot the saving gains ###############################
# obtain the optimal cost dataframe
optimal_cost_eff_df <- optimal_power_df_eff |>
  dplyr::group_by(assay, cost_ratio, primer_efficiency_bound) |>
  dplyr::filter(total_cost == min(total_cost)) |>
  dplyr::ungroup() |>
  dplyr::mutate(minimum_fold_change = as.numeric(as.character(minimum_fold_change))) |> 
  dplyr::select(cost_ratio, total_cost, assay, primer_efficiency_bound) |>
  dplyr::mutate(assay = ifelse(assay == "Perturb-seq", "Perturb_seq", "TAP_seq")) |>
  tidyr::pivot_wider(names_from = c("assay"), values_from = c("total_cost")) |> 
  dplyr::mutate(Perturb_TAP_cost_ratio = Perturb_seq / TAP_seq,
                primer_filtering = primer_efficiency_bound) |>
  dplyr::left_join(parameter_grid |> dplyr::rename(cost_ratio = ratio), 
                   by = c("cost_ratio", "primer_efficiency_bound"))

# Get data point at cost_ratio = 0.23 for annotation
annotation_point <- optimal_cost_eff_df |>
  dplyr::filter(arm_ratio, round(cost_ratio,2) == 0.23) |>
  dplyr::mutate(label = "10x/Illumina\nprices 2025")

cost_over_cost_ratio_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_ratio) |>
  ggplot(aes(x = 1 / cost_ratio, y = 1 - 1 / Perturb_TAP_cost_ratio)) +
  geom_point(color = "black") +
  geom_line(color = "black", linewidth = 1) +
  geom_text_repel(
    data = annotation_point,
    aes(label = label),
    size = 4,
    nudge_x = -0.4,
    nudge_y = 0.15,
    direction = "both",
    box.padding = 0.5,
    point.padding = 0.5,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  # <-- new annotation
  annotate(
    "text",
    x = Inf, y = 0.25,                    # near bottom of auto-scaled range
    label = "Primer efficiency > 40%",
    hjust = 1.1, vjust = 0,              # text sits above the anchor
    colour = "grey40",
    size = 4
  ) +
  scale_x_log10() +
  labs(
    x = "Sequencing cost / Library prep cost",
    y = "TAP-seq cost savings",
    title = "Cost saving versus relative cost"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5, size = 16),
    plot.tag     = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.9)
  )

# Cost over primer efficiency plot
cost_over_primer_efficiency_plot <- optimal_cost_eff_df |>
  dplyr::filter(arm_primer_efficiency_bound) |>
  ggplot(aes(x = primer_filtering, y = 1 - 1 / Perturb_TAP_cost_ratio)) +
  geom_point(color = "black") +
  geom_line(color = "black", linewidth = 1) +
  labs(x = "Minimum primer efficiency", y = "TAP-seq cost savings",
       title = "Cost saving versus primer efficiency", tag = "e") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.9)
  )

# add inset plot
source("reproduction-code-prep/case-study/TAP-seq-case-study/saturation_relative_expression_comparison.R")
cost_over_primer_efficiency_plot <- cost_over_primer_efficiency_plot + 
  patchwork::inset_element(
    relative_expression_profile_inset +
      labs(tag = NULL) +                      
      theme(plot.tag = element_blank()),      
    left = 0.4, bottom = 0, right = 1, top = 0.58,
    align_to = "panel"
  ) 

# save the combined plot 
saveRDS(cost_over_cost_ratio_plot, sprintf("%s/cell_read_cost_cost_ratio_plotting.rds", dir_to_results))
saveRDS(cost_over_primer_efficiency_plot, sprintf("%s/cell_read_cost_primer_eff_plotting.rds", dir_to_results))