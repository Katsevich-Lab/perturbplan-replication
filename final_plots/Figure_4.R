# This is a script creating Figure 3 in the manuscript
library(tibble)      # tribble()
library(ggplot2)     # ggplotify depends on it anyway
library(gridExtra)   # tableGrob()
library(patchwork)   # plot_layout()
library(cowplot)     # get_plot_component()
library(scales)      # trans_breaks

# specify the results folder
dir_to_results <- "reproduction-code-prep/case-study/perturb-seq-case-study/results/"

# define the plots_dir
plots_dir <- "reproduction-code-prep/final_plots/figures"
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# specify power target
power_target <- 0.8

################### Create Figure 3(a), a table ################################
# ── data ─────────────────────────────────────────────────────────────────────
tbl_df <- tribble(
  ~Metric,                    ~Morris, ~Gasperini, ~Replogle,
  "Number of targets",        "600",   "6K",       "10K",
  "MOI",                      "13",    "28",       "1",
  "Number of cells",          "46K",   "200K",     "2M",
  "Reads per cell",           "60K",   "34K",      "15K",
  "Cells per target",         "1K",    "1k",      "200"
)

# ── build the table grob ─────────────────────────────────────────────────────
ggplot_colors <- scales::hue_pal()(3)  # Get default ggplot2 colors for 3 categories

# ── build the table grob ─────────────────────────────────────────────────────
tbl_grob <- tableGrob(
  tbl_df, rows = NULL,
  theme = ttheme_minimal(
    base_size = 22,
    core    = list(fg_params = list(hjust = .5, x = .5)),    # centre the body
    colhead = list(                                           # colored header bar
      fg_params = list(col = "white", fontface = "bold"),
      bg_params = list(fill = c("white", ggplot_colors), col = NA)
    )
  )
)

############################## Create Figure 4(a) ##############################
# load the power over TPM 
power_over_tpm <- readRDS(paste0(dir_to_results, "Morris/power_over_tpm_df.rds")) |> 
  dplyr::mutate(assay = "Morris") |>
  dplyr::bind_rows(readRDS(paste0(dir_to_results, "Gasperini/power_over_tpm_df.rds")) |>
                     dplyr::mutate(assay = "Gasperini")) |>
  dplyr::bind_rows(readRDS(paste0(dir_to_results, "Replogle/power_over_tpm_df.rds")) |>
                     dplyr::mutate(assay = "Replogle")) |> 
  dplyr::mutate(assay = factor(assay, levels = c("Morris", "Gasperini", "Replogle"))) 

# do the plotting
power_over_tpm_plot <- power_over_tpm |> 
  dplyr::filter(TPM_threshold <= 1000) |>
  ggplot(aes(x = TPM_threshold, y = overall_power, color = assay)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1) +
  scale_x_log10() + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x = "TPM analysis threshold", y = "Overall power", title = sprintf("Varying TPM with FC = %d%%", round((1 - unique(power_over_tpm$minimum_fold_change)) * 100))) +
  theme_classic() +
  scale_color_manual(values = c("Morris" = "#fdae61", "Gasperini" = "#abd9e9", "Replogle" = "#4575b4")) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x  = element_text(size = 14)  # ← strip label size
  )

############################## Create Figure 4(b) ##############################
# load the power over FC
power_over_fc <- readRDS(paste0(dir_to_results, "Morris/power_over_fc_df.rds")) |> 
  dplyr::mutate(assay = "Morris") |>
  dplyr::bind_rows(readRDS(paste0(dir_to_results, "Gasperini/power_over_fc_df.rds")) |>
                     dplyr::mutate(assay = "Gasperini")) |>
  dplyr::bind_rows(readRDS(paste0(dir_to_results, "Replogle/power_over_fc_df.rds")) |>
                     dplyr::mutate(assay = "Replogle")) |> 
  dplyr::mutate(assay = factor(assay, levels = c("Morris", "Gasperini", "Replogle")))

# do the plotting
power_over_fc_plot <- power_over_fc |> 
  dplyr::mutate(minimum_fold_change = round((1 - minimum_fold_change) * 100, 1)) |> 
  dplyr::filter(minimum_fold_change <= 50) |>
  ggplot(aes(x = minimum_fold_change, y = overall_power, color = assay)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x = "Minimum fold change (percent)", y = "Overall power", title = "Varying FC with all expressed genes") +
  theme_classic() +
  scale_color_manual(values = c("Morris" = "#fdae61", "Gasperini" = "#abd9e9", "Replogle" = "#4575b4")) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x  = element_text(size = 14)  # ← strip label size
  )

# obtain the shared legend
shared_legend <- get_plot_component(power_over_fc_plot, "guide-box", return_all = TRUE)[[3]]

########################## combine plot to get Figure 2 (a,b) ##################
combined_power_over_tpm_FC <- ((power_over_tpm_plot / power_over_fc_plot) & 
                                 theme(legend.position = "none")) + 
  plot_annotation(
    title = "Overall power vs. TPM and Minimum FC",
    theme = theme(plot.title = element_text(hjust = 0.5, size  = 16))
  )

############################# Figure 4 (c, d) #####################################
opt_power_tpm <- readRDS(paste0(dir_to_results, "Morris/power_opt_tpm_df.rds")) |>
  dplyr::mutate(assay = "Morris") |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Gasperini/power_opt_tpm_df.rds")) |>
      dplyr::mutate(assay = "Gasperini")
  ) |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Replogle/power_opt_tpm_df.rds")) |>
      dplyr::mutate(assay = "Replogle")
  ) |> 
  dplyr::mutate(assay = factor(assay, levels = c("Morris", "Gasperini", "Replogle")))

opt_cost_tpm <- readRDS(paste0(dir_to_results, "Morris/cost_opt_tpm_df.rds")) |>
  dplyr::mutate(assay = "Morris") |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Gasperini/cost_opt_tpm_df.rds")) |>
      dplyr::mutate(assay = "Gasperini")
  ) |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Replogle/cost_opt_tpm_df.rds")) |>
      dplyr::mutate(assay = "Replogle")
  ) |> 
  dplyr::mutate(assay = factor(assay, levels = c("Morris", "Gasperini", "Replogle")))

# obtain optimal cost df
tpm_to_plot <- c(10, 50, 200)
optimal_cost_tpm_df <- opt_power_tpm |> dplyr::filter(total_cost == minimum_cost) 

# plot the reads over effect size
cell_read_plot_varying_tpm <- optimal_cost_tpm_df |>
  dplyr::mutate(TPM_threshold = factor(TPM_threshold)) |>
  dplyr::filter(TPM_threshold %in% as.character(tpm_to_plot)) |>
  ggplot(aes(x = num_captured_cells, y = sequenced_reads_per_cell, color = assay)) +
  scale_y_continuous(limits = c(30000, 135000), 
                     labels = c("40K", "60K", "80K", "100K", "120K"), 
                     breaks = c(40000, 60000, 80000, 100000, 120000)) +
  scale_x_log10(limits = c(9000, 10000000), 
                labels = c("10K", "100K", "1M", "10M"), 
                breaks = c(10000, 100000, 1000000, 10000000)) +
  geom_point(size = 1) +
  geom_line(size = 1) +
  geom_point(
    data = optimal_cost_tpm_df |> dplyr::filter(TPM_threshold %in% as.character(tpm_to_plot)),
    aes(x = num_captured_cells,     
        y = sequenced_reads_per_cell,            
        shape = TPM_threshold),
    size  = 3
  ) +
  guides(colour  = "none", linetype = guide_legend(order = 1, override.aes = list(colour = "black"))) +
  labs(x = "Number of total cells", y = "Reads per cell", title = "Optimal (cells, reads) vs. TPM",
       shape = "TPM threshold") +
  theme_classic() +
  guides(shape = guide_legend(
    title.position = "left",   # title on same row, at left
    label.position = "right",
    nrow = 1,                  # all keys in one row
    byrow = TRUE
  )) +
  scale_color_manual(values = c("Morris" = "#fdae61", "Gasperini" = "#abd9e9", "Replogle" = "#4575b4")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = c(0.55, 0.9),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x  = element_text(size = 14)  # ← strip label size
  )

# Cost over TPM plot
cost_tpm_plot <- optimal_cost_tpm_df |>
  dplyr::mutate(TPM_threshold = as.numeric(as.character(TPM_threshold))) |>
  ggplot(aes(x = TPM_threshold, y = minimum_cost, color = assay)) +
  geom_point() +
  scale_y_log10(
    breaks = 10^(3:6),                                   # 10^3, 10^4, 10^5, 10^6
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL                                   # optional: hide minor ticks
  ) +
  scale_x_log10() +
  geom_line(size = 1) +
  labs(x = "TPM analysis threshold", y = "Optimal cost ($)", title = "Optimal cost vs. TPM threshold") +
  theme_classic() +
  scale_color_manual(values = c("Morris" = "#fdae61", "Gasperini" = "#abd9e9", "Replogle" = "#4575b4")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x  = element_text(size = 14)  # ← strip label size
  )

############################# Figure 4(e, f) #####################################
opt_power_fc <- readRDS(paste0(dir_to_results, "Morris/power_opt_fc_df.rds")) |>
  dplyr::mutate(assay = "Morris") |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Gasperini/power_opt_fc_df.rds")) |>
      dplyr::mutate(assay = "Gasperini")
  ) |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Replogle/power_opt_fc_df.rds")) |>
      dplyr::mutate(assay = "Replogle")
  ) |> 
  dplyr::mutate(assay = factor(assay, levels = c("Morris", "Gasperini", "Replogle")))

opt_cost_fc <- readRDS(paste0(dir_to_results, "Morris/cost_opt_fc_df.rds")) |>
  dplyr::mutate(assay = "Morris") |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Gasperini/cost_opt_fc_df.rds")) |>
      dplyr::mutate(assay = "Gasperini")
  ) |>
  dplyr::bind_rows(
    readRDS(paste0(dir_to_results, "Replogle/cost_opt_fc_df.rds")) |>
      dplyr::mutate(assay = "Replogle")
  ) |>
  dplyr::mutate(assay = factor(assay, levels = c("Morris", "Gasperini", "Replogle")))

# obtain optimal cost df
fc_to_plot <- c(20, 30, 40)
optimal_cost_fc_df <- opt_power_fc |> dplyr::filter(total_cost == minimum_cost)
optimal_cost_fc_df |> 
  dplyr::filter(minimum_fold_change == 20) |> 
  dplyr::select(num_captured_cells, total_cost)

# plot the reads over effect size
df_plot <- optimal_cost_fc_df |>
  dplyr::mutate(
    minimum_fold_change = as.character(minimum_fold_change)
  ) |>
  dplyr::filter(minimum_fold_change %in% as.character(fc_to_plot)) |>
  dplyr::mutate(
    minimum_fold_change = factor(
      minimum_fold_change,
      levels = as.character(sort(as.numeric(fc_to_plot)))
    )
  )

cell_read_plot_varying_fc <- ggplot(df_plot, aes(x = num_captured_cells, y = sequenced_reads_per_cell, color = assay)) +
  geom_point() +
  geom_line(size = 1) +
  geom_point(
    aes(shape = minimum_fold_change),
    size = 3
  ) +
  scale_y_continuous(
    limits = c(40000, 130000),
    labels = c("40K", "60K", "80K", "100K", "120K"),
    breaks = c(40000, 60000, 80000, 100000, 120000)
  ) +
  scale_x_log10(
    limits = c(5000, 20000000),
    labels = c("10K", "100K", "1M", "10M"),
    breaks = c(10000, 100000, 1000000, 10000000)
  ) +
  labs(
    x = "Number of total cells", y = "Reads per cell",
    shape = "Minimum FC",
    title = "Optimal (cells, reads) vs. minimum FC"
  ) +
  theme_classic() +
  scale_color_manual(values = c("Morris" = "#fdae61", "Gasperini" = "#abd9e9", "Replogle" = "#4575b4")) +
  guides(colour = "none") +
  guides(shape = guide_legend(
    title.position = "left",
    label.position = "right",
    nrow = 1,
    byrow = TRUE
  )) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = c(0.55, 0.1),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x = element_text(size = 14)
  )

# Cost over FC
cost_fc_plot <- optimal_cost_fc_df |>
  dplyr::mutate(minimum_fold_change = as.numeric(as.character(minimum_fold_change))) |>
  ggplot(aes(x = minimum_fold_change, y = minimum_cost, color = assay)) +
  geom_point() +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),      # 1, 10, 100, …
    labels = trans_format("log10", math_format(10^.x))      # 10^0, 10^1, …
  ) +
  geom_line(size = 1) +
  labs(x = "Minimum fold change (percent)", y = "Optimal cost ($)", title = "Optimal cost vs. minimum FC") +
  theme_classic() +
  scale_color_manual(values = c("Morris" = "#fdae61", "Gasperini" = "#abd9e9", "Replogle" = "#4575b4")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text.x  = element_text(size = 14)  # ← strip label size
  )

######################## Collate the final plots ###############################
# combine the plot for optimal (cells, reads)
combined_cell_read_plot <- cell_read_plot_varying_tpm / cell_read_plot_varying_fc
combined_cost_plot <- cost_tpm_plot / cost_fc_plot

# wrap the plots
Figure_4 <- wrap_plots(
  combined_power_over_tpm_FC,
  plot_spacer(),                # <- empty space
  combined_cell_read_plot,
  plot_spacer(),                # <- empty space
  combined_cost_plot,
  widths = c(1, 0.3, 1, 0.1, 1),
  ncol = 5, byrow = FALSE
) 

# save the plot
ggsave(file.path(plots_dir, "figure_4.pdf"), Figure_4, width = 14, height = 8)
