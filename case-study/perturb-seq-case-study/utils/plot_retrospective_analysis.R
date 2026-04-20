library(ggplot2)
library(patchwork)

# do the power over effect size plot
varying_fc_plot <- power_over_fc_df |> 
  dplyr::mutate(minimum_fold_change = round((1 - minimum_fold_change) * 100, 1)) |> 
  ggplot(aes(x = minimum_fold_change, y = overall_power)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              se = FALSE, color = "blue", size = 1) +
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = "Minimum effect size (percent)", y = "Overall power", title = "Varying fold change") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  )

# do the power over effect size plot
varying_tpm_plot <- power_over_tpm_df |> 
  ggplot(aes(x = TPM_threshold, y = overall_power)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              se = FALSE, color = "blue", size = 1) +
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = "TPM analysis threshold", y = "Overall power", title = "Varying TPM threshold") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  )

# plot the power curve
varying_fc_tpm_plot <- varying_fc_tpm_df |> 
  ggplot(aes(x = minimum_fold_change, y = TPM_threshold)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              se = FALSE, color = "blue", size = 1) +
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = "Minimum effect size (percent)", y = "TPM analysis threshold") +
  theme_bw() +
  labs(title = sprintf("Target power at %.2f", power_target)) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  )

############################# combine and save plots ###########################
merged_plot <- varying_fc_plot + varying_tpm_plot + varying_fc_tpm_plot
