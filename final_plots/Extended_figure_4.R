# This is a script creating supplementary pilot figure
library(tibble)      # tribble()
library(ggplot2)     # ggplotify depends on it anyway
library(gridExtra)   # tableGrob()
library(ggplotify)   # as.ggplot()
library(patchwork)   # plot_layout()
library(cowplot)     # get_plot_component()
library(scales)      # trans_breaks
library(dplyr)
library(perturbplan)

# specify the results folder
dir_to_results <- "reproduction-code-prep/case-study/perturb-seq-case-study/results/"

# define the plots_dir
plots_dir <- "reproduction-code-prep/final_plots/figures"
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# specify power target
power_target <- 0.8

# specify colormapping
color_map <- c(
  "K562 (10x)" = "#2166ac",
  "K562 (Gasperini)" = "#b2182b"
)

############################## Create Figure supp (a) ##############################
# load the power over TPM
power_over_tpm <- readRDS(paste0(dir_to_results, "Gasperini/power_over_tpm_df.rds")) |>
  dplyr::mutate(assay = "K562 (Gasperini)") |>
  dplyr::bind_rows(readRDS(paste0(dir_to_results, "K562_10x/power_over_tpm_df.rds")) |>
                     dplyr::mutate(assay = "K562 (10x)")) |>
  dplyr::mutate(assay = factor(assay, levels = c("K562 (Gasperini)", "K562 (10x)")))

# do the plotting
power_over_tpm_plot <- power_over_tpm |>
  dplyr::filter(TPM_threshold <= 1000) |>
  dplyr::mutate(type = "Varying TPM with FC = 25%") |>
  ggplot(aes(x = TPM_threshold, y = overall_power, color = assay)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_color_manual(values = color_map) + 
  facet_wrap(.~type) +
  labs(x = "TPM analysis threshold", y = "Overall power") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x  = element_text(size = 15)  # ← strip label size
  )

############################## Create Figure supp (b) ##############################
# load the power over FC
power_over_fc <- readRDS(paste0(dir_to_results, "Gasperini/power_over_fc_df.rds")) |>
  dplyr::mutate(assay = "K562 (Gasperini)") |>
  dplyr::bind_rows(readRDS(paste0(dir_to_results, "K562_10x/power_over_fc_df.rds")) |>
                     dplyr::mutate(assay = "K562 (10x)")) |>
  dplyr::mutate(assay = factor(assay, levels = c("K562 (Gasperini)", "K562 (10x)")))

# do the plotting
power_over_fc_plot <- power_over_fc |>
  dplyr::mutate(minimum_fold_change = round((1 - minimum_fold_change) * 100, 1)) |>
  dplyr::filter(minimum_fold_change <= 50) |>
  dplyr::mutate(type = "Varying FC with all expressed genes") |>
  ggplot(aes(x = minimum_fold_change, y = overall_power, color = assay)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1) +
  facet_wrap(.~type) +
  scale_color_manual(values = color_map) + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x = "Minimum fold change (percent)", y = "Overall power") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x  = element_text(size = 15)  # ← strip label size
  )

# obtain the shared legend
shared_legend <- get_plot_component(power_over_fc_plot, "guide-box", return_all = TRUE)[[3]]

########################## combine plot to get Figure supp (a,b) ##################
combined_power_over_tpm_FC <- ((power_over_tpm_plot + power_over_fc_plot) &
                                 theme(legend.position = "none")) +
  plot_annotation(
    title = "Overall power vs. TPM and Minimum FC",
    theme = theme(plot.title = element_text(hjust = 0.5, size  = 16))
  )

############################# Figure supp (c, d) #####################################
# Only create these plots if the required files exist for both datasets
if(file.exists(paste0(dir_to_results, "Gasperini/power_opt_tpm_df.rds")) &&
   file.exists(paste0(dir_to_results, "K562_10x/power_opt_tpm_df.rds"))) {

  opt_power_tpm <- readRDS(paste0(dir_to_results, "Gasperini/power_opt_tpm_df.rds")) |>
    dplyr::mutate(assay = "K562 (Gasperini)") |>
    dplyr::bind_rows(
      readRDS(paste0(dir_to_results, "K562_10x/power_opt_tpm_df.rds")) |>
        dplyr::mutate(assay = "K562 (10x)")
    ) |>
    dplyr::mutate(assay = factor(assay, levels = c("K562 (Gasperini)", "K562 (10x)")))

  opt_cost_tpm <- readRDS(paste0(dir_to_results, "Gasperini/cost_opt_tpm_df.rds")) |>
    dplyr::mutate(assay = "K562 (Gasperini)") |>
    dplyr::bind_rows(
      readRDS(paste0(dir_to_results, "K562_10x/cost_opt_tpm_df.rds")) |>
        dplyr::mutate(assay = "K562 (10x)")
    ) |>
    dplyr::mutate(assay = factor(assay, levels = c("K562 (Gasperini)", "K562 (10x)")))

  # obtain optimal cost df
  tpm_to_plot <- c(10, 50, 200)
  optimal_cost_tpm_df <- opt_power_tpm |> dplyr::filter(total_cost == minimum_cost)

  # plot the reads over effect size
  cell_read_plot_varying_tpm <- optimal_cost_tpm_df |>
    dplyr::mutate(TPM_threshold = factor(TPM_threshold)) |>
    dplyr::filter(TPM_threshold %in% as.character(tpm_to_plot)) |>
    dplyr::mutate(type = "Optimal (cells, reads) vs. TPM") |>
    ggplot(aes(x = num_captured_cells, y = sequenced_reads_per_cell, color = assay)) +
    facet_wrap(.~type) +
    scale_color_manual(values = color_map) + 
    scale_y_continuous(limits = c(30000, 130000), labels = comma_format()) +
    scale_x_log10(labels = comma_format()) +
    geom_point(size = 1) +
    geom_line(size = 1) +
    geom_point(
      data = optimal_cost_tpm_df |> dplyr::filter(TPM_threshold %in% as.character(tpm_to_plot)),
      aes(x = num_captured_cells,
          y = sequenced_reads_per_cell,
          shape = TPM_threshold),
      size  = 3
    ) +
    labs(x = "Number of total cells", y = "Reads per cell",
         shape = "TPM threshold") +
    guides(colour = "none", linetype = guide_legend(order = 1, override.aes = list(colour = "black"))) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.position = c(0.84, 0.25),  # Bottom right position
      plot.title = element_text(hjust = 0.5, size = 16),
      strip.text.x  = element_text(size = 15)  # ← strip label size
    )

  # Cost over TPM plot
  cost_tpm_plot <- optimal_cost_tpm_df |>
    dplyr::mutate(TPM_threshold = as.numeric(as.character(TPM_threshold))) |>
    dplyr::mutate(type = "Optimal cost vs. TPM threshold") |>
    ggplot(aes(x = TPM_threshold, y = minimum_cost, color = assay)) +
    facet_wrap(.~type) +
    geom_point() +
    scale_y_log10(
      breaks = 10^(3:6),                                   # 10^3, 10^4, 10^5, 10^6
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      minor_breaks = NULL                                   # optional: hide minor ticks
    ) +
    scale_x_log10() +
    scale_color_manual(values = color_map) + 
    geom_line(size = 1) +
    labs(x = "TPM analysis threshold", y = "Optimal cost ($)") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 16),
      strip.text.x  = element_text(size = 15)  # ← strip label size
    )
} else {
  # Create empty placeholders if files don't exist
  cell_read_plot_varying_tpm <- ggplot() +
    theme_void() +
    ggtitle("Data not available for optimal (cells, reads) vs. TPM")
  cost_tpm_plot <- ggplot() +
    theme_void() +
    ggtitle("Data not available for optimal cost vs. TPM")
}

############################# Figure supp (e, f) #####################################
# Only create these plots if the required files exist for both datasets
if(file.exists(paste0(dir_to_results, "Gasperini/power_opt_fc_df.rds")) &&
   file.exists(paste0(dir_to_results, "K562_10x/power_opt_fc_df.rds"))) {

  opt_power_fc <- readRDS(paste0(dir_to_results, "Gasperini/power_opt_fc_df.rds")) |>
    dplyr::mutate(assay = "K562 (Gasperini)") |>
    dplyr::bind_rows(
      readRDS(paste0(dir_to_results, "K562_10x/power_opt_fc_df.rds")) |>
        dplyr::mutate(assay = "K562 (10x)")
    ) |>
    dplyr::mutate(assay = factor(assay, levels = c("K562 (Gasperini)", "K562 (10x)")))

  opt_cost_fc <- readRDS(paste0(dir_to_results, "Gasperini/cost_opt_fc_df.rds")) |>
    dplyr::mutate(assay = "K562 (Gasperini)") |>
    dplyr::bind_rows(
      readRDS(paste0(dir_to_results, "K562_10x/cost_opt_fc_df.rds")) |>
        dplyr::mutate(assay = "K562 (10x)")
    ) |>
    dplyr::mutate(assay = factor(assay, levels = c("K562 (Gasperini)", "K562 (10x)")))

  # obtain optimal cost df
  fc_to_plot <- c(20, 30, 40)
  optimal_cost_fc_df <- opt_power_fc |> dplyr::filter(total_cost == minimum_cost)

  # plot the reads over effect size
  cell_read_plot_varying_fc <- optimal_cost_fc_df |>
    dplyr::filter(minimum_fold_change %in% fc_to_plot) |>
    dplyr::mutate(minimum_fold_change = as.numeric(as.character(minimum_fold_change))) |>
    dplyr::mutate(type = "Optimal (cells, reads) vs. minimum FC") |>
    ggplot(aes(x = num_captured_cells, y = sequenced_reads_per_cell, color = assay)) +
    facet_wrap(.~type) +
    geom_point() +
    scale_y_continuous(limits = c(30000, 130000), labels = comma_format()) +
    geom_line(size = 1) +
    scale_color_manual(values = color_map) + 
    scale_x_log10(labels = comma_format()) +
    geom_point(
      data = optimal_cost_fc_df |> dplyr::filter(minimum_fold_change %in% as.character(fc_to_plot)),
      aes(x = num_captured_cells, y = sequenced_reads_per_cell, shape = minimum_fold_change),
      size  = 3
    ) +
    labs(x = "Number of total cells", y = "Reads per cell", shape = "Minimum fold change") +
    guides(colour = "none", linetype = guide_legend(order = 1, override.aes = list(colour = "black"))) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.position = c(0.7, 0.3),  # Bottom right position, moved a bit left
      plot.title = element_text(hjust = 0.5, size = 16),
      strip.text.x  = element_text(size = 15)
    )

  # Cost over FC
  cost_fc_plot <- optimal_cost_fc_df |>
    dplyr::mutate(minimum_fold_change = as.numeric(as.character(minimum_fold_change))) |>
    dplyr::mutate(type = "Optimal cost vs. minimum FC") |>
    ggplot(aes(x = minimum_fold_change, y = minimum_cost, color = assay)) +
    facet_wrap(.~type) +
    geom_point() +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),      # 1, 10, 100, …
      labels = trans_format("log10", math_format(10^.x))      # 10^0, 10^1, …
    ) +
    scale_color_manual(values = color_map) + 
    geom_line(size = 1) +
    labs(x = "Minimum fold change (percent)", y = "Optimal cost ($)") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 16),
      strip.text.x  = element_text(size = 15)  # ← strip label size
    )
} else {
  # Create empty placeholders if files don't exist
  cell_read_plot_varying_fc <- ggplot() +
    theme_void() +
    ggtitle("Data not available for optimal (cells, reads) vs. FC")
  cost_fc_plot <- ggplot() +
    theme_void() +
    ggtitle("Data not available for optimal cost vs. FC")
}

######################## Collate the final plots ###############################
# Remove legends from power and cost plots, but keep internal legends in cell_read plots
power_over_tpm_plot_no_legend <- power_over_tpm_plot + theme(legend.position = "none")
power_over_fc_plot_no_legend <- power_over_fc_plot + theme(legend.position = "none")
cell_read_plot_varying_tpm_no_legend <- cell_read_plot_varying_tpm  # Keep internal shape legend
cost_tpm_plot_no_legend <- cost_tpm_plot + theme(legend.position = "none")
cell_read_plot_varying_fc_no_legend <- cell_read_plot_varying_fc  # Keep internal shape legend
cost_fc_plot_no_legend <- cost_fc_plot + theme(legend.position = "none")

# Create empty plot for spacing
empty_plot <- ggplot() + theme_void()

########################### Create top 2 plots #################################
data("K562_Gasperini")
data("K562_10x")

# ============================================================
# Panel a/b: case study (K562 10x) TPM + dispersion
# ============================================================
# default expression stats from package (K562)
default_expr <- perturbplan:::get_pilot_data_from_package("K562")$baseline_expression_stats %>%
  dplyr::filter(!is.na(relative_expression)) %>%
  rename(default_expr = relative_expression, default_size = expression_size)

# custom expression
custom_expr <- K562_10x$baseline_expression_stats %>%
  dplyr::filter(!is.na(relative_expression)) %>%
  rename(custom_expr = relative_expression, custom_size = expression_size)

panel_a <- default_expr %>%
  inner_join(custom_expr, by = "response_id") %>%
  mutate(default_expr = default_expr * 1e6, custom_expr = custom_expr * 1e6) %>%
  ggplot(aes(default_expr, custom_expr)) +
  geom_point(alpha = 0.1) +
  geom_abline() +
  scale_x_log10(breaks = c(1, 10, 30, 100, 1000)) +
  scale_y_log10(breaks = c(1, 10, 30, 100, 1000)) +
  labs(x = "K562 (Gasperini)", y = "K562 (10x)", title = "Gene-level TPM comparison") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

panel_b <- default_expr %>%
  inner_join(custom_expr, by = "response_id") %>%
  ggplot(aes(default_size, custom_size)) +
  geom_point(alpha = 0.1) +
  geom_abline() +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),      # 1, 10, 100, …
    labels = trans_format("log10", math_format(10^.x))      # 10^0, 10^1, …
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),      # 1, 10, 100, …
    labels = trans_format("log10", math_format(10^.x))      # 10^0, 10^1, …
  ) +
  labs(x = "K562 (Gasperini)", y = "K562 (10x)", title = "Disperison (size) comparison") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# Create the main plot grid: 2 columns x 4 rows
final_plots <- wrap_plots(
  panel_a,                                      # a (row 1, col 1) - relative expression
  panel_b,                                      # b (row 1, col 2) - size parameter
  power_over_tpm_plot_no_legend,                # c (row 1, col 1) - power vs TPM
  power_over_fc_plot_no_legend,                 # d (row 1, col 2) - power vs FC
  cell_read_plot_varying_tpm_no_legend,         # e (row 2, col 1) - (cell, reads) vs TPM
  cell_read_plot_varying_fc_no_legend,          # f (row 2, col 2) - (cell, reads) vs FC
  cost_tpm_plot_no_legend,                      # g (row 3, col 1) - cost vs TPM
  cost_fc_plot_no_legend,                       # h (row 3, col 2) - cost vs FC
  ncol = 2,
  nrow = 4,
  byrow = TRUE  # byrow = TRUE to fill by row, labels go a-f across rows
) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag          = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.9),         # inside panel, top-left
    plot.title        = element_text(size = 18, hjust = 0.5, vjust = 1)
  )

# Create bottom row with legend in the middle
bottom_row <- cowplot::plot_grid(
  NULL,
  shared_legend,
  NULL,
  ncol = 3,
  rel_widths = c(1, 1, 1)
)

# Combine main plots with legend row using cowplot
Figure_supp_pilot <- cowplot::plot_grid(
  final_plots,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 0.04)
)

################################################################################
########################### Generate Figure S3 #################################
################################################################################
# Save as extended_figure_4.pdf
ggsave(file.path(plots_dir, "extended_figure_4.pdf"), Figure_supp_pilot, width = 10, height = 13.5)
