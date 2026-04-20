# This is a script comparing TAP-seq and Perturb-seq
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)
library(patchwork)

# source necessary files
source("reproduction-code-prep/case-study/helper.R")

# load TAP-seq pilot data
pilot_data <- readRDS("reproduction-code-prep/case-study/auxiliary-files/discovery_es_at_scale.rds")

# set the primer efficiency threshold: lower bound on share-normalized
# TAP/Perturb ratio. 0.4 corresponds to "40% primer efficiency".
primer_efficiency_threshold <- 0.4

# set the TPM threshold
UMI_threshold <- 1

# extract the data (no primer-efficiency filtering at this stage)
extracted_pilot_data <- get_baseline_expression(primer_efficiency_threshold = 0, pilot_data = pilot_data)

# path to power results
dir_to_results <- "reproduction-code-prep/case-study/TAP-seq-case-study/results/"
if(!dir.exists(dir_to_results)){
  dir.create(dir_to_results, recursive = TRUE)
}

############################## Obtain fixed parameters #########################
fixed_parameters <- readRDS("reproduction-code-prep/case-study/fixed_parameters.rds")
mapping_efficiency_perturb <- 1
mapping_efficiency_tap <- fixed_parameters$mapping_efficiency_tap / fixed_parameters$mapping_efficiency

# calculate representation of target genes (TAP-seq genes) in perturb-seq K562 data
k562_data <- perturbplan:::get_pilot_data_from_package("K562")
tap_genes <- extracted_pilot_data$baseline_expression_stats_tap$response_id
representation_target_genes <- k562_data$baseline_expression_stats |>
  dplyr::filter(response_id %in% tap_genes) |>
  dplyr::pull(relative_expression) |>
  sum()
mapping_efficiency_perturb_adjusted <- mapping_efficiency_perturb * representation_target_genes

######################## compare sequencing saturation #########################
# construct plotting dataframe
reads_per_cell <- 10^{seq(1, 6, length.out = 100)}
read_umi_df <- data.frame(
  reads_per_cell_sequenced = reads_per_cell,
  reads_per_cell_mapped = reads_per_cell * mapping_efficiency_tap,
  assay = "TAP-seq"
) |>
  dplyr::mutate(
    library_size = perturbplan:::fit_read_UMI_curve_cpp(
      reads_per_cell = reads_per_cell_mapped,
      rSAC_fn_wrapper = extracted_pilot_data$library_parameters_tap
    ),
    umis_per_cell_goi = library_size
  ) |>
  dplyr::bind_rows(
    data.frame(
      reads_per_cell_sequenced = reads_per_cell,
      reads_per_cell_mapped = reads_per_cell * mapping_efficiency_perturb_adjusted,
      assay = "Perturb-seq"
    ) |>
      dplyr::mutate(
        library_size = perturbplan:::fit_read_UMI_curve_cpp(
          reads_per_cell = reads_per_cell_sequenced * mapping_efficiency_perturb,
          rSAC_fn_wrapper = k562_data$library_parameters
        ),
        umis_per_cell_goi = library_size * representation_target_genes
  )
)

# compute reads/cell to reach 80 percent saturation
saturation_reads_98_percent <- read_umi_df |> 
  dplyr::group_by(assay) |> 
  dplyr::mutate(saturation_98 = max(0.98 * max(library_size))) |> 
  dplyr::mutate(library_size_diff = abs(library_size - saturation_98)) |> 
  dplyr::filter(library_size_diff == min(library_size_diff)) |> 
  dplyr::ungroup() 

# do the plot
saturation_comparison <- read_umi_df |>
  ggplot(aes(x = reads_per_cell_sequenced, y = umis_per_cell_goi, color = assay)) +
  geom_line(size = 1) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),      # 1, 10, 100, …
    labels = trans_format("log10", math_format(10^.x))      # 10^0, 10^1, …
  ) +
  scale_y_continuous() +
  ## show only linetype legend inside plot, hide color legend
  guides(colour = "none", linetype = guide_legend(order = 1)) +
  labs(x = "Sequenced reads per cell", y = "UMIs per cell on target",
       color = "Assay",
       title = "Saturation curve comparison",
       tag = "a") +
  theme_classic() +
  scale_color_manual(values = c("Perturb-seq" = "#a50026", "TAP-seq" = "#313695")) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.9)
  )

# create barplot inset for mapping efficiencies
efficiency_df <- data.frame(
  assay = c("TAP-seq", "Perturb-seq"),
  efficiency = c(mapping_efficiency_tap, mapping_efficiency_perturb_adjusted)
)

# get the colors used in the main plot
assay_colors <- scales::hue_pal()(2)  # default ggplot2 colors
names(assay_colors) <- c("Perturb-seq", "TAP-seq")

efficiency_barplot <- ggplot(efficiency_df, aes(x = assay, y = efficiency, fill = assay)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", efficiency * 100)),
            vjust = -0.5, size = 4, color = "black") +
  scale_fill_manual(values = assay_colors) +
  labs(title = "Fraction of reads \non targeted genes") +
  theme_classic() +
  scale_fill_manual(values = c("Perturb-seq" = "#a50026", "TAP-seq" = "#313695")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.tag = element_blank()
  ) +
expand_limits(y = 1.2 * max(efficiency_df$efficiency)) 

# combine main plot with inset
saturation_comparison <- saturation_comparison +
  patchwork::inset_element(
    efficiency_barplot +
      labs(tag = NULL) +                      
      theme(plot.tag = element_blank()),      
    left = 0.02, bottom = 0.3, right = 0.4, top = 0.9,
    align_to = "panel"
  ) +
  patchwork::plot_annotation(tag_levels = NULL)

########################## Compare relative expression #########################
# transform the relative expression to TPM
representation_target_genes <- sum(extracted_pilot_data$baseline_expression_stats_perturb$relative_expression)

# construct relative expression distribution
relative_expression_df <- extracted_pilot_data$baseline_expression_stats_tap |> 
  dplyr::rename(TAP_seq = UMI_per_cell) |> 
  dplyr::left_join(
    extracted_pilot_data$baseline_expression_stats_perturb |> 
      dplyr::select(UMI_per_cell, response_id) |>
      dplyr:::rename(Perturb_seq = UMI_per_cell) |>
      dplyr::distinct(), 
    by = "response_id"
  )  

# Compute the boundary slope in (TAP_UMI, Perturb_UMI) space corresponding to
# the share-normalized TAP/Perturb threshold:
#   e_j = (TAP_UMI_j / tap_total) / (Perturb_UMI_j / perturb_total) >= thr
# <=> Perturb_UMI_j <= (perturb_total / (thr * tap_total)) * TAP_UMI_j
tap_total      <- sum(relative_expression_df$TAP_seq)
perturb_total  <- sum(relative_expression_df$Perturb_seq, na.rm = TRUE)
boundary_slope <- perturb_total / (primer_efficiency_threshold * tap_total)

# construct the label for the linear function annotation
labels_df <- tibble::tribble(
  ~x,     ~y,        ~label,
  300,  0.1,  sprintf("%d%% primer \nefficiency", round(primer_efficiency_threshold * 100)),
)

# plot the truncation scatter plot. The reference line is the share-normalized
# "100% primer efficiency" locus: Perturb_UMI = (perturb_total / tap_total) * TAP_UMI.
relative_expression_profile <- relative_expression_df |>
  ggplot(aes(TAP_seq, Perturb_seq)) +
  geom_point() +
  geom_abline(slope = 1, intercept = log10(perturb_total / tap_total)) +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  labs(x = "TAP-seq", y = "Perturb-seq",
       title = "UMIs/cell at saturation comparison") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    strip.text   = element_text(size = 14),
    plot.title   = element_text(hjust = 0.5, size = 16),
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.95)
  )

# create the inset plot
relative_expression_profile_inset <- relative_expression_df |>
  mutate(in_band = !is.na(Perturb_seq) & Perturb_seq <= boundary_slope * TAP_seq) |>
  ggplot(aes(TAP_seq, Perturb_seq, colour = in_band)) +
  geom_point() +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"), guide = "none") +
  geom_abline(slope = 1, intercept = log10(perturb_total / tap_total)) +
  ## share-normalized threshold line
  annotate(
    "segment",
    x    = 0.01,
    xend = 100,
    y    = boundary_slope * 0.01,
    yend = boundary_slope * 100,
    colour   = "red",
    linetype = "dashed"
  ) +
  ## ——— add the curve labels ———
  geom_text(data = labels_df,
            aes(x, y, label = label),
            angle   = 0,        # aligns with slope-1 lines on log–log axes
            hjust   = 1,
            size    = 3,
            color = "red",
            inherit.aes = FALSE) +
  scale_x_log10(
    breaks = c(0.01, 1, 100),
    labels = c("0.01", "1", "100")
  ) +
  scale_y_log10(
    breaks = c(0.01, 1, 100),
    labels = c("0.01", "1", "100")
  ) +
  labs(x = "TAP-seq", y = "Perturb-seq", title = "UMIs/cell at saturation") +
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    plot.title   = element_text(hjust = 0.5, size = 12),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.tag = element_blank()
  ) 

########################### Save the plot separately ###########################
# combine the plots (no additional legend at bottom - legend will be added in Figure4.R)
saveRDS(saturation_comparison, paste0(dir_to_results, "saturation_plotting.rds"))
saveRDS(relative_expression_profile, paste0(dir_to_results, "UMI_comparison_plotting.rds"))
