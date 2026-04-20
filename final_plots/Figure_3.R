# This is a script attempting collate the benchmarking and validation results for Figure 2 in the paper
library(ggplot2)
library(patchwork)
library(cowplot)
library(stringr)

# define the plots_dir
plots_dir <- "reproduction-code-prep/final_plots/figures"
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# define the top folder for results
results_dir <- "reproduction-code-prep"

# load the plotting object for benchmarking (Figure 2(a))
simulation_benchmark <- readRDS(sprintf("%s/simulation-benchmarking/benchmarking_plot.rds", results_dir))

# load the plotting object for read-UMI curve (Figure 2(b))
read_umi_validation <- readRDS(sprintf("%s/saturation-curve-fitting/results/saturation_curve_plot.rds", results_dir))

# load the plotting object for Gasperini (high effect size) (Figure 2(c))
Gasperini_high_effect_size <- readRDS(sprintf("%s/realdata-validation-pipeline/PerturbPlan/results/Gasperini_gene_analysis_plot.rds", results_dir))
Gasperini_legend <- readRDS(sprintf("%s/realdata-validation-pipeline/PerturbPlan/results/Gasperini_gene_analysis_legend_plot.rds", results_dir))

# load the plotting object for Ray (low effect size) (Figure 2(d))
Ray_low_effect_size <- readRDS(sprintf("%s/realdata-validation-pipeline/PerturbPlan/results/Ray_enhancer_analysis_plot.rds", results_dir))

####################### collate the results #####################################
# 1. Build the 2 × 2 block **once**, tell patchwork to tag just this level
four_panels <- (
  (simulation_benchmark | read_umi_validation) /
    ((Gasperini_high_effect_size + theme(legend.position = "none")) |
       (Ray_low_effect_size        + theme(legend.position = "none")))
) +
  plot_layout(heights   = c(1, 3.5)) +     # <-- stop tagging deeper
  plot_annotation(                     # automatic a, b, c, d
    tag_levels = "a"
  ) & 
  theme(
    strip.background = element_blank(),
    plot.tag          = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.9)      # top-left inside each panel
  )

# 2. Add the legend row underneath (legend is *outside* the tagged block)
Figure_3 <- cowplot::plot_grid(
  four_panels,
  Gasperini_legend,            # extracted earlier with get_legend()
  ncol        = 1,
  rel_heights = c(1, 0.03)
)

# save the plot
ggsave(file.path(plots_dir, "figure_3.pdf"), Figure_3, width = 15, height = 15)
