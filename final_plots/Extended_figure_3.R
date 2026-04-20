library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(cowplot)
library(minpack.lm)
library(preseqR)
source("reproduction-code-prep/saturation-curve-fitting/helper.R")

# -----------------------------
# Output location
# -----------------------------
plots_dir <- "reproduction-code-prep/final_plots/figures"
out_pdf   <- file.path(plots_dir, "extended_figure_3.pdf")
num_SRR_in_use <- 1

# specify the results folder
dir_to_results <- "reproduction-code-prep/saturation-curve-fitting/results"

# add plotting theme
common_validation_theme <- function() {
  theme_bw() +
    theme(
      axis.title.x = element_text(size = 12),
      axis.text.x  = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.y  = element_text(size = 10),
      legend.text  = element_text(size = 12),
      legend.title = element_blank(),
      strip.text   = element_text(size = 12),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA)
    )
}

# ============================================================
# Part A: Validation plots (T_CD8 + K562)  -> panels a/b
# ============================================================
# load data directory/
plot_data_validation <- readRDS(sprintf("%s/saturation_curve_comparison.rds", dir_to_results))
plot_data_validation <- plot_data_validation |>
  dplyr::mutate(
    source = dplyr::case_when(
      source %in% c("Perturbplan fit", "PerturbPlan fit", "PerturbPlan (scPower) fit") ~ "PerturbPlan (preseqR) fit",
      TRUE ~ source
    )
  )

# add lintype
linetype_map_validation <- c(
  "PerturbPlan (preseqR) fit" = "dashed",
  "Downsampled real data" = "solid",
  "scPower fit"           = "solid"
)

# add color
color_map_validation <- c(
  "PerturbPlan (preseqR) fit" = "#67001f",
  "Downsampled real data" = "#f4a582",
  "scPower fit"           = "#4393c3"
)

# K562 validation plot
k562_validation_plot <- plot_data_validation %>%
  dplyr::filter(cell_type == "K562 (Gasperini)") %>%
  ggplot(aes(x = reads_per_cell, y = library_size, color = source, linetype = source)) +
  geom_line(linewidth = 1.25) +
  scale_linetype_manual(values = linetype_map_validation) +
  scale_color_manual(values = color_map_validation) +
  facet_wrap(~ cell_type, scales = "free") +
  labs(x = "Mapped reads per cell", y = "UMIs per cell", color = NULL, linetype = NULL) +
  expand_limits(y = 0) +
  common_validation_theme()

# TCD8 validation plot
tcd8_validation_plot <- plot_data_validation %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "T CD8 (Shifrut)", "CD8+ T (Shifrut)", cell_type)) |>
  dplyr::filter(cell_type == "CD8+ T (Shifrut)") %>%
  ggplot(aes(x = reads_per_cell, y = library_size, color = source, linetype = source)) +
  geom_line(linewidth = 1.25) +
  scale_linetype_manual(values = linetype_map_validation) +
  scale_color_manual(values = color_map_validation) +
  facet_wrap(~ cell_type, scales = "free") +
  labs(x = "Mapped reads per cell", y = "UMIs per cell", color = NULL, linetype = NULL) +
  expand_limits(y = 0) +
  common_validation_theme()

# shared legend from the K562 plot
shared_legend <- get_plot_component(k562_validation_plot, "guide-box", return_all = TRUE)[[3]]
k562_no_legend <- k562_validation_plot + theme(legend.position = "none")
tcd8_no_legend <- tcd8_validation_plot + theme(legend.position = "none")

# combine plots
top_row <- cowplot::plot_grid(
  k562_no_legend, tcd8_no_legend,
  ncol = 2,
  labels = c("a", "b"),
  label_size = 12,
  label_fontface = "bold"
)
top_row_with_legend <- cowplot::plot_grid(
  top_row, shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

# ============================================================
# Part B: Varying sequencing saturation plot (K562 varying_sd) -> panel c
# ============================================================
# load data file
read_umi_df_vsd <- readRDS(sprintf("%s/read_umi_df_varying_sd.rds", dir_to_results))
read_umi_grid_df_vsd <- readRDS(sprintf("%s/saturation_curve_comparison_varying_sequencing_depth.rds", dir_to_results))
library_parameters_vsd <- readRDS(sprintf("%s/library_parameters_on_1_SRRs_varying_sd.rds", dir_to_results))
hline_y <- library_parameters_vsd[[length(library_parameters_vsd)]]$UMI_per_cell_at_saturation

# Get unique downsampling ratios from the data
unique_ds_ratios <- unique(read_umi_grid_df_vsd$downsampling_ratio)
unique_ds_ratios <- sort(unique_ds_ratios)

# color map for panels c and d (unified color scheme) - now using downsampling_ratio as key
manual_saturation_colors <- c(
  "0.01" = "#b2182b",
  "0.03" = "#d6604d",
  "0.1" = "#f4a582",
  "0.3" = "#92c5de",
  "1" = "#2166ac"
)

# linetype map - now using downsampling_ratio as key
linetype_map_vsd <- c(
  "0.01" = "solid",
  "0.03" = "solid",
  "0.1" = "solid",
  "0.3" = "solid",
  "1" = "dashed"
)

# Calculate x position for full dataset vertical line
full_reads_per_cell <- read_umi_df_vsd %>%
  dplyr::filter(downsampling_ratio == 1) %>%
  dplyr::pull(reads_per_cell) %>%
  mean()

# do plots over different sequencing depth
panel_c <- ggplot(read_umi_grid_df_vsd, aes(
      x = reads_per_cell / hline_y,
      y = library_size / hline_y,
      color = factor(round(downsampling_ratio, 3)),
      linetype = factor(round(downsampling_ratio, 3))
    )) +
  geom_line(
    alpha = 0.75,
    linewidth = 1.25
  ) +
  # Add vertical line marking full dataset position
  geom_vline(xintercept = full_reads_per_cell / hline_y, color = "gray40", linetype = "dotted", linewidth = 0.8) +
  scale_linetype_manual(values = linetype_map_vsd, guide = "none") +
  scale_x_log10() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(
    values = manual_saturation_colors,
    labels = function(x) {
      sapply(x, function(ds_ratio) {
        ds_ratio <- as.numeric(ds_ratio)
        if (abs(ds_ratio - 1) < 0.001) {
          return("100%")
        } else if (abs(ds_ratio - 0.3) < 0.001) {
          return("30%")
        } else if (abs(ds_ratio - 0.1) < 0.001) {
          return("10%")
        } else if (abs(ds_ratio - 0.03) < 0.001) {
          return("3%")
        } else if (abs(ds_ratio - 0.01) < 0.001) {
          return("1%")
        } else {
          return(sprintf("%.1f%%", ds_ratio * 100))
        }
      })
    }
  ) +
  geom_hline(yintercept = 1, color = "blue", linetype = "dashed") +
  annotate(
    "text",
    x = 0.1, y = 1,
    label  = sprintf("UMIs at saturation \nin K562 (Gasperini) = %.0f", hline_y),
    hjust  = 0.8, vjust = 1.2, size = 3, colour = "blue"
  ) +
  labs(
    x = "Reads per cell / UMIs at saturation",
    y = "Library size / UMIs at saturation",
    color = "Reference dataset size (fraction of full dataset)",
    title = "Curves learned on downsampled K562 (Gasperini)"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 12, hjust = 0.5, vjust = 2),
    axis.title.x    = element_text(size = 12),
    axis.text.x     = element_text(size = 10),
    axis.title.y    = element_text(size = 12),
    axis.text.y     = element_text(size = 10),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = "bottom"
  )

# ============================================================
# Part C: Varying sequencing saturation plot (TCD8 varying_sd) -> panel d
# ============================================================
# load data file
read_umi_df_vsd_tcd8 <- readRDS(sprintf("%s/read_umi_df_varying_sd_tcd8.rds", dir_to_results))
read_umi_grid_df_vsd_tcd8 <- readRDS(sprintf("%s/saturation_curve_comparison_varying_sequencing_depth_tcd8.rds", dir_to_results))
library_parameters_vsd_tcd8 <- readRDS(sprintf("%s/library_parameters_tcd8_varying_sd.rds", dir_to_results))
hline_y_tcd8 <- library_parameters_vsd_tcd8[[length(library_parameters_vsd_tcd8)]]$UMI_per_cell_at_saturation

# Get unique downsampling ratios from the data
unique_ds_ratios_tcd8 <- unique(read_umi_grid_df_vsd_tcd8$downsampling_ratio)
unique_ds_ratios_tcd8 <- sort(unique_ds_ratios_tcd8)

# Calculate x position for full dataset vertical line
full_reads_per_cell_tcd8 <- read_umi_df_vsd_tcd8 %>%
  dplyr::filter(downsampling_ratio == 1) %>%
  dplyr::pull(reads_per_cell) %>%
  mean()

# do plots over different sequencing depth
panel_d <- ggplot(read_umi_grid_df_vsd_tcd8, aes(
      x = reads_per_cell / hline_y_tcd8,
      y = library_size / hline_y_tcd8,
      color = factor(round(downsampling_ratio, 3)),
      linetype = factor(round(downsampling_ratio, 3))
    )) +
  geom_line(
    alpha = 0.75,
    linewidth = 1.25
  ) +
  # Add vertical line marking full dataset position
  geom_vline(xintercept = full_reads_per_cell_tcd8 / hline_y_tcd8, color = "gray40", linetype = "dotted", linewidth = 0.8) +
  scale_linetype_manual(values = linetype_map_vsd, guide = "none") +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(
    values = manual_saturation_colors,
    labels = function(x) {
      sapply(x, function(ds_ratio) {
        ds_ratio <- as.numeric(ds_ratio)
        if (abs(ds_ratio - 1) < 0.001) {
          return("100%")
        } else if (abs(ds_ratio - 0.3) < 0.001) {
          return("30%")
        } else if (abs(ds_ratio - 0.1) < 0.001) {
          return("10%")
        } else if (abs(ds_ratio - 0.03) < 0.001) {
          return("3%")
        } else if (abs(ds_ratio - 0.01) < 0.001) {
          return("1%")
        } else {
          return(sprintf("%.1f%%", ds_ratio * 100))
        }
      })
    }
  ) +
  geom_hline(yintercept = 1, color = "purple", linetype = "dashed") +
  annotate(
    "text",
    x = 0.1, y = 1,
    label  = sprintf("UMIs at saturation \nin CD8+ T (Shifrut) = %.0f", hline_y_tcd8),
    hjust  = 0.6, vjust = 1.2, size = 3, colour = "purple"
  ) +
  labs(
    x = "Reads per cell / UMIs at saturation",
    y = "Library size / UMIs at saturation",
    color = "Reference dataset size (fraction of full dataset)",
    title = "Curves learned on downsampled CD8+ T (Shifrut)"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 12, hjust = 0.5, vjust = 2),
    axis.title.x    = element_text(size = 12),
    axis.text.x     = element_text(size = 10),
    axis.title.y    = element_text(size = 12),
    axis.text.y     = element_text(size = 10),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = "bottom"
  )

# ============================================================
# Part E: Extrapolation error analysis for K562 -> panel e
# ============================================================
# Load true data for K562 (this provides the x-grid)
read_umi_df_k562_true <- readRDS(sprintf("%s/read_umi_df_k562.rds", dir_to_results))

# Get unique downsampling ratios
ds_ratios_k562 <- sort(unique(read_umi_df_vsd$downsampling_ratio))

# For each downsampling ratio, generate predictions on the true data's x-grid
extrapolation_error_k562 <- dplyr::bind_rows(
  lapply(seq_along(ds_ratios_k562), function(i) {
    ds_ratio <- ds_ratios_k562[i]

    # Get library parameters for this downsampling ratio
    lib_params <- library_parameters_vsd[[i]]

    # Get the training data for this downsampling ratio (to find initial point)
    training_data <- read_umi_df_vsd %>%
      dplyr::filter(downsampling_ratio == ds_ratio)

    # Get the maximum reads_per_cell for this downsampling ratio (initial training point)
    initial_reads <- max(training_data$reads_per_cell, na.rm = TRUE)

    # Generate predictions on the true data's x-grid using the library parameters
    predicted_library_size <- perturbplan:::fit_read_UMI_curve_cpp(
      read_umi_df_k562_true$reads_per_cell,
      lib_params
    )

    # Create result data frame
    result <- data.frame(
      reads_per_cell = read_umi_df_k562_true$reads_per_cell,
      actual_library_size = read_umi_df_k562_true$library_size,
      predicted_library_size = predicted_library_size,
      downsampling_ratio = ds_ratio,
      initial_reads = initial_reads
    ) %>%
      dplyr::mutate(
        extrapolation_factor = reads_per_cell / initial_reads,
        error_pct = (predicted_library_size - actual_library_size) / actual_library_size * 100
      )

    return(result)
  })
)

panel_e <- extrapolation_error_k562 %>%
  dplyr::filter(extrapolation_factor >= 1, round(downsampling_ratio, 3) != 1) %>%
  ggplot(aes(x = extrapolation_factor, y = error_pct,
             color = factor(round(downsampling_ratio, 3)))) +
  geom_line(linewidth = 1.25) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(
    values = manual_saturation_colors,
    labels = function(x) {
      sapply(x, function(ds_ratio) {
        ds_ratio <- as.numeric(ds_ratio)
        if (abs(ds_ratio - 1) < 0.001) {
          return("100%")
        } else if (abs(ds_ratio - 0.3) < 0.001) {
          return("30%")
        } else if (abs(ds_ratio - 0.1) < 0.001) {
          return("10%")
        } else if (abs(ds_ratio - 0.03) < 0.001) {
          return("3%")
        } else if (abs(ds_ratio - 0.01) < 0.001) {
          return("1%")
        } else {
          return(sprintf("%.1f%%", ds_ratio * 100))
        }
      })
    }
  ) +
  labs(
    x = "Extrapolation factor",
    y = "Extrapolation error (%)",
    color = "Reference dataset size (fraction of full dataset)",
    title = "Extrapolation error for K562 (Gasperini)"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 12, hjust = 0.5, vjust = 2),
    axis.title.x    = element_text(size = 12),
    axis.text.x     = element_text(size = 10),
    axis.title.y    = element_text(size = 12),
    axis.text.y     = element_text(size = 10),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = "bottom"
  )

# ============================================================
# Part F: Extrapolation error analysis for CD8+ T -> panel f
# ============================================================
# Load true data for CD8+ T (this provides the x-grid)
read_umi_df_tcd8_true <- readRDS(sprintf("%s/read_umi_df_tcd8.rds", dir_to_results))

# Get unique downsampling ratios
ds_ratios_tcd8 <- sort(unique(read_umi_df_vsd_tcd8$downsampling_ratio))

# For each downsampling ratio, generate predictions on the true data's x-grid
extrapolation_error_tcd8 <- dplyr::bind_rows(
  lapply(seq_along(ds_ratios_tcd8), function(i) {
    ds_ratio <- ds_ratios_tcd8[i]

    # Get library parameters for this downsampling ratio
    lib_params <- library_parameters_vsd_tcd8[[i]]

    # Get the training data for this downsampling ratio (to find initial point)
    training_data <- read_umi_df_vsd_tcd8 %>%
      dplyr::filter(downsampling_ratio == ds_ratio)

    # Get the maximum reads_per_cell for this downsampling ratio (initial training point)
    initial_reads <- max(training_data$reads_per_cell, na.rm = TRUE)

    # Generate predictions on the true data's x-grid using the library parameters
    predicted_library_size <- perturbplan:::fit_read_UMI_curve_cpp(
      read_umi_df_tcd8_true$reads_per_cell,
      lib_params
    )

    # Create result data frame
    result <- data.frame(
      reads_per_cell = read_umi_df_tcd8_true$reads_per_cell,
      actual_library_size = read_umi_df_tcd8_true$library_size,
      predicted_library_size = predicted_library_size,
      downsampling_ratio = ds_ratio,
      initial_reads = initial_reads
    ) %>%
      dplyr::mutate(
        extrapolation_factor = reads_per_cell / initial_reads,
        error_pct = (predicted_library_size - actual_library_size) / actual_library_size * 100
      )

    return(result)
  })
)

panel_f <- extrapolation_error_tcd8 %>%
  dplyr::filter(extrapolation_factor >= 1, round(downsampling_ratio, 3) != 1) %>%
  ggplot(aes(x = extrapolation_factor, y = error_pct,
             color = factor(round(downsampling_ratio, 3)))) +
  geom_line(linewidth = 1.25) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(
    values = manual_saturation_colors,
    labels = function(x) {
      sapply(x, function(ds_ratio) {
        ds_ratio <- as.numeric(ds_ratio)
        if (abs(ds_ratio - 1) < 0.001) {
          return("100%")
        } else if (abs(ds_ratio - 0.3) < 0.001) {
          return("30%")
        } else if (abs(ds_ratio - 0.03) < 0.001) {
          return("3%")
        } else if (abs(ds_ratio - 0.1) < 0.001) {
          return("10%")
        } else if (abs(ds_ratio - 0.01) < 0.001) {
          return("1%")
        } else {
          return(sprintf("%.1f%%", ds_ratio * 100))
        }
      })
    }
  ) +
  labs(
    x = "Extrapolation factor",
    y = "Extrapolation error (%)",
    color = "Reference dataset size (fraction of full dataset)",
    title = "Extrapolation error for CD8+ T (Shifrut)"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 12, hjust = 0.5, vjust = 2),
    axis.title.x    = element_text(size = 12),
    axis.text.x     = element_text(size = 10),
    axis.title.y    = element_text(size = 12),
    axis.text.y     = element_text(size = 10),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = "bottom"
  )

# ============================================================
# Combine and save (ONLY output)
# ============================================================
# Create 3x2 layout: K562 (left column) vs CD8+ T (right column)
# Row 1: Original data fit (panels a and b)
# Row 2: Downsampled data comparison (panels c and d)
# Row 3: Extrapolation error (panels e and f)

# Extract legends
# Legend for a & b (validation plots)
legend_validation <- get_plot_component(k562_validation_plot, "guide-box", return_all = TRUE)[[3]]

# Legend for c, d, e, f (pilot dataset size)
legend_downsampling <- get_plot_component(panel_c, "guide-box", return_all = TRUE)[[3]]

# Remove legends from individual panels
panel_c_no_legend <- panel_c + theme(legend.position = "none")
panel_d_no_legend <- panel_d + theme(legend.position = "none")
panel_e_no_legend <- panel_e + theme(legend.position = "none")
panel_f_no_legend <- panel_f + theme(legend.position = "none")

# Create row 1 with validation plots
row_1 <- cowplot::plot_grid(
  k562_no_legend, tcd8_no_legend,
  ncol = 2,
  labels = c("a", "b"),
  label_size = 16,
  label_fontface = "bold"
)

row_1_with_legend <- cowplot::plot_grid(
  row_1,
  legend_validation,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

# Create rows 2 and 3 with downsampling plots
rows_2_3 <- cowplot::plot_grid(
  panel_c_no_legend, panel_d_no_legend,
  panel_e_no_legend, panel_f_no_legend,
  ncol = 2,
  nrow = 2,
  labels = c("c", "d", "e", "f"),
  label_size = 16,
  label_fontface = "bold",
  rel_heights = c(1, 0.8)
)

rows_2_3_with_legend <- cowplot::plot_grid(
  rows_2_3,
  legend_downsampling,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

# Combine all rows
extended_figure_3_combined <- cowplot::plot_grid(
  row_1_with_legend,
  rows_2_3_with_legend,
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 2)
)

# Save the plot
save_plot(extended_figure_3_combined, out_pdf, width = 9, height = 11)
