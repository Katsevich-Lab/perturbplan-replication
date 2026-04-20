# load path to raw data directory
library(perturbplan)
library(ggplot2)
library(scales)

# load learned pilot data
data("K562_Gasperini")

# make directory
dir_to_results <- "reproduction-code-prep/saturation-curve-fitting/results"
if(!dir.exists(dir_to_results)){
  dir.create(dir_to_results)
}

# read data 
read_umi_learning <- readRDS(sprintf("%s/read_umi_learning_10x_Gasperini.rds", dir_to_results))
read_umi_downsampling <- readRDS(sprintf("%s/read_umi_downsampling_10x_Gasperini.rds", dir_to_results))

# define auxiliary vairables
hline_y        <- K562_Gasperini$library_parameters$UMI_per_cell_at_saturation
max_reads      <- max(read_umi_learning$reads_per_cell)
end_of_curve_y <- max(read_umi_learning |> 
                        dplyr::filter(reads_per_cell == max_reads) |> 
                        dplyr::pull(library_size))

# do line plot
saturation_curve_plot <- read_umi_downsampling |> 
  dplyr::mutate(assay = ifelse(assay == "10x example", "10x", "Gasperini")) |>
  ggplot(aes(x = reads_per_cell, y = library_size, color = assay)) +
  geom_line(size = 1.25) +
  scale_y_continuous(labels = c("20k", "40k", "60k"), breaks = c(20000, 40000, 60000)) +
  scale_x_continuous(labels = comma_format()) +
  geom_line(
    data        = read_umi_learning,
    mapping     = aes(x = reads_per_cell, y = library_size, linetype = fitted_approach),
    color = "grey50",
    alpha       = 0.5,      # ← transparency here
    size        = 1.25,       # ← line thickness (default is 0.5 mm)
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = hline_y, color = "blue", linetype = "dashed") +
  
  ## ── NEW: annotate the two reference geoms ─────────────────────────────
  annotate(
    "text",
    x      = max_reads * 0.80,   # 80 % of the way across
    y      = hline_y,
    label  = "Estimated total UMIs/cell",
    hjust  = 0.6, vjust = 1.2,    # tweak so it sits just above the line
    size   = 4, colour = "blue"
  ) +
  annotate(
    "text",
    x      = max_reads * 0.80,
    y      = end_of_curve_y,
    label  = "Learned saturation curve",
    hjust  = 0.75, vjust = -0.6,
    size   = 4
  ) +
  ## ──────────────────────────────────────────────────────────────────────
  
  labs(x = "Reads per cell",
       y = "UMIs per cell",
       title = "Saturation curve fitting", 
       color = "Data",
       linetype = "Learned on") +
  theme_classic() +
  scale_color_manual(values = c("10x" = "#a50026", "Gasperini" = "#313695")) +
  theme(
    legend.position = c(0.65, 0.14),
    legend.box = "horizontal",
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    axis.title.x   = element_text(size = 14),
    axis.text.x    = element_text(size = 12),
    axis.title.y   = element_text(size = 14),
    axis.text.y    = element_text(size = 12),
    legend.spacing.x = unit(-0.02, "cm")
  )

# save the plot object
saveRDS(saturation_curve_plot, sprintf("%s/saturation_curve_plot.rds", dir_to_results))
