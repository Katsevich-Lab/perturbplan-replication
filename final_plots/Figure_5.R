# This is a script collating the subplots for Figure 4 in the main text
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(perturbplan)
library(ggplotify)

# load the plotting object
dir_to_rds <- "reproduction-code-prep/case-study/TAP-seq-case-study/results/"

# define the plots_dir
plots_dir <- "reproduction-code-prep/final_plots/figures"
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# load the subplot for figure 2(a)
saturation_comparison <- readRDS(paste0(dir_to_rds, "saturation_plotting.rds"))

# load the subplot for figure 2(b, c)
power_versus_tpm <- readRDS(paste0(dir_to_rds, "power_versus_TPM_plotting.rds"))
power_versus_fc <- readRDS(paste0(dir_to_rds, "power_versus_FC_plotting.rds"))

# load the subplot for figure 2(d)
UMI_comparison <- readRDS(paste0(dir_to_rds, "UMI_comparison_plotting.rds"))

# load the subplot for figure 2(e, f)
TAP_Perturb_saving_plot_cost_ratio <- readRDS(paste0(dir_to_rds, "cell_read_cost_cost_ratio_plotting.rds"))
TAP_Perturb_saving_plot_primer_eff <- readRDS(paste0(dir_to_rds, "cell_read_cost_primer_eff_plotting.rds"))

# load the shared legend
shared_legend <- readRDS(paste0(dir_to_rds, "shared_legend_plotting.rds"))

###################### combined plots ##########################################
three_panels <- (((saturation_comparison + labs(tag = "a")) / (UMI_comparison + labs(tag = "b"))) | 
                   ((power_versus_tpm + labs(tag = "c")) / (power_versus_fc + labs(tag = "d"))) |
                   ((TAP_Perturb_saving_plot_primer_eff + labs(tag = "e")) / (TAP_Perturb_saving_plot_cost_ratio + labs(tag = "f")))
)  +
  plot_annotation(tag_levels = NULL) 

# add the shared legend
Figure_5 <- cowplot::plot_grid(
  three_panels,
  shared_legend,            # extracted earlier with get_legend()
  ncol        = 1,
  rel_heights = c(1, 0.06)
)

# save the plot
ggsave(file.path(plots_dir, "figure_5.pdf"), Figure_5, width = 13, height = 8)
