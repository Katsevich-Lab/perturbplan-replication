# This is a script plotting the results for TAP-seq versus Perturb-seq in retrospective analyses
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)  # for unit()

# do the power over effect size plot
varying_fc_plot <- power_over_fc_df |>
  dplyr::filter(minimum_fold_change >= 0.75) |>
  dplyr::mutate(
    minimum_fold_change = round((1 - minimum_fold_change) * 100, 1)
  ) |>
  dplyr::filter(minimum_fold_change <= 50) |>
  ggplot(aes(x = minimum_fold_change, y = overall_power, color = assay)) +
  # <-- new annotation
  annotate(
    "text",
    x = Inf, y = 0.05,                    # corner of the panel
    label = sprintf("UMIs/cell > %d", UMI_threshold_default),
    hjust = 1.1, vjust = 1.5,            # nudge inward a bit
    colour = "grey40",
    size = 4
  ) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 1) +
  geom_hline(yintercept = power_target, color = "grey40", linetype = "dashed") +
  labs(x = "Minimum effect size (percent)", y = "Overall power", 
       title = "Varying minimum effect size", tag = "d") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  scale_color_manual(values = c("Perturb-seq" = "#a50026", "TAP-seq" = "#313695")) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.9)
  )

# do the power over effect size plot
varying_tpm_plot <- power_over_tpm_df |>
  ggplot(aes(x = UMI_threshold, y = overall_power, color = assay)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = FALSE, linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "TAP-seq expression threshold (UMIs/cell)", y = "Overall power", 
       title = "Varying TAP-seq UMIs/cell", tag = "c") +
  theme_classic() +
  scale_color_manual(values = c("Perturb-seq" = "#a50026", "TAP-seq" = "#313695")) +
  geom_hline(yintercept = power_target, color = "grey40", linetype = "dashed") +
  annotate(
    "text",
    x      = 0.5,   # 80 % of the way across
    y      = power_target,
    label  = sprintf("Power target = %.1f", power_target),
    hjust  = 0.6, vjust = 1.2,    # tweak so it sits just above the line
    size   = 4, colour = "grey40"
  ) +
  annotate(
    "text",
    x = Inf, y = 0.05,                    # corner of the panel
    label = sprintf("Effect size = %d%%", round((1 - minimum_fold_change) * 100)),
    hjust = 1.1, vjust = 1.5,            # nudge inward a bit
    colour = "grey40",
    size = 4
  ) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.9)
  )

# save the Power over Fold Change and TPM analysis threshold
saveRDS(varying_fc_plot & theme(legend.position = "none"), paste0(dir_to_results, "power_versus_FC_plotting.rds"))
saveRDS(varying_tpm_plot & theme(legend.position = "none"), paste0(dir_to_results, "power_versus_TPM_plotting.rds"))

# obtain the legend and save it
shared_legend <- get_plot_component(varying_fc_plot + theme(legend.position = "bottom", legend.direction = "horizontal"), "guide-box", return_all = TRUE)[[3]]
saveRDS(shared_legend, paste0(dir_to_results, "shared_legend_plotting.rds"))
