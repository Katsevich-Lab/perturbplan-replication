# This is a Rscript checkinf the monotonicity of FDP curve
library(dplyr) 
library(tidyr) 
library(ggplot2)
library(patchwork)
library(stringr)
library(cowplot)
library(scales)

# Assign arguments to variables
experiment <- "Gasperini"
positive_proportion <- 0.01
alpha <- 0.01
window_def <- "gene_analysis"
guides_per_target <- 2
downsampled_guides_per_target <- 2
cutoff_grid <- exp(seq(log(1e-5), log(1e-1), length.out = 500))

# source the necessary files
source("~/.Rprofile")
source("reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R")
source("reproduction-code-prep/realdata-validation-pipeline/run_fc_distr_construction.R")

# specify intermediate files
intermediate_files_folder <- "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/"

# specify the fc_tail_list (replacing fc_spread_list)
fc_tail_list <- c("light_tail", "heavy_tail")

########################## 1. Compute library size #############################
# fit S-M curve
source("reproduction-code-prep/realdata-validation-pipeline/run_fitting_S_M_curve.R")

# load the parameter grid
plan_parameter_grid <- readRDS(sprintf("%s/%s/real_data_parameter_grid.rds", intermediate_files_folder, experiment))

# obtain the total UMIs by fitting a nonlinear model
library_size_estimate <- plan_parameter_grid |>
  dplyr::rename(num_total_cells = N, reads_per_cell = R) |>
  dplyr::select(num_total_cells, reads_per_cell, grid_id)

############################### 2. Loop over (R, N) pairs ######################
# Initialize both result dataframes
analytic_results <- NULL
for (grid_id_cur in library_size_estimate$grid_id) {
  
  # subset the dataframe according to grid_id
  reads_cells_df <- library_size_estimate |> dplyr::filter(grid_id == grid_id_cur)
  
  # extract the number of cells, reads_per_cell and library parameters
  num_total_cells <- reads_cells_df$num_total_cells
  reads_per_cell <- reads_cells_df$reads_per_cell 
  library_parameters <- fitted_S_M_curve
  
  ############################### 3. Compute power for both approaches ###########
  # Initialize power dataframes 
  overall_power_max <- data.frame(matrix(NA, ncol = length(fc_tail_list)))
  overall_power_median <- data.frame(matrix(NA, ncol = length(fc_tail_list)))
  colnames(overall_power_median) <- fc_tail_list
  colnames(overall_power_max) <- fc_tail_list
  
  # Initialize storage for num_trt_cells by tail type
  num_trt_cells_by_tail <- setNames(rep(NA, length(fc_tail_list)), fc_tail_list)
  
  for (fc_tail in fc_tail_list) {
    
    # Load unified dataframe with both cell count methods
    discovery_set <- readRDS(paste0(intermediate_files_folder, experiment, sprintf("/fc_distr_%s_%s.rds", fc_tail, window_def)))
    
    ############################### 3.1 Compute common quantities (same for both approaches) #####################
    # Compute fold_change and baseline expression (identical for both approaches)
    fold_change_summary <- discovery_set |> dplyr::summarise(max_fc = max(fold_change), median_fc = median(fold_change))
    fold_change_sd <- unique(discovery_set$gRNA_sd)
    
    # Compute library size (same for both approaches)
    library_size <- perturbplan::fit_read_UMI_curve_cpp(
      reads_per_cell = reads_per_cell * targeted_genes_reads_ratio,
      rSAC_fn_wrapper = library_parameters
    )
    
    # Compute baseline gene expression (same for both approaches)
    custom_baseline_data <- list(baseline_expression = 
                                   list(baseline_expression = discovery_set |> 
                                          dplyr::select(response_id, relative_expression, expression_size)))
    fc_expression_df_max <- perturbplan::extract_fc_expression_info(minimum_fold_change = fold_change_summary$max_fc,
                                                                    gRNA_variability = fold_change_sd,
                                                                    B = 1000, TPM_threshold = 0,
                                                                    custom_pilot_data = custom_baseline_data,
                                                                    gRNAs_per_target = downsampled_guides_per_target)
    fc_expression_df_median <- perturbplan::extract_fc_expression_info(minimum_fold_change = fold_change_summary$median_fc,
                                                                       gRNA_variability = fold_change_sd,
                                                                       B = 1000, TPM_threshold = 0,
                                                                       custom_pilot_data = custom_baseline_data,
                                                                       gRNAs_per_target = downsampled_guides_per_target)
    
    ############################### 3.2 Compute approach-specific cell counts #####################
    # Regular approach cell counts
    num_trt_cells <- discovery_set |>
      dplyr::select(grna_target, num_trt_cells_full, num_total_plan_cells) |>
      dplyr::distinct() |>
      dplyr::mutate(num_cells = round(num_trt_cells_full * num_total_cells / num_total_plan_cells)) |> 
      dplyr::summarise(num_trt_cells = round(mean(num_cells) * downsampled_guides_per_target / guides_per_target)) |>
      dplyr::pull()
    num_cntrl_cells <- num_total_cells - num_trt_cells
    
    # Store num_trt_cells for this tail type
    num_trt_cells_by_tail[fc_tail] <- num_trt_cells
    
    ############################ 3.3 Compute means and sds #####################
    # Regular approach
    max_mean_sd_list <- perturbplan:::compute_monte_carlo_teststat_cpp(
      fc_expression_df = fc_expression_df_max$fc_expression_df,
      library_size = library_size, 
      num_trt_cells = num_trt_cells, 
      num_cntrl_cells = num_cntrl_cells
    )
    median_mean_sd_list <- perturbplan:::compute_monte_carlo_teststat_cpp(
      fc_expression_df = fc_expression_df_median$fc_expression_df,
      library_size = library_size, 
      num_trt_cells = num_trt_cells, 
      num_cntrl_cells = num_cntrl_cells
    )

    # Append results
    analytic_results <- rbind(analytic_results,
                              # Lower bound for this tail
                              data.frame(
                                reads_per_cell = unname(reads_per_cell),
                                num_total_cells = unname(num_total_cells), 
                                num_trt_cells = unname(num_trt_cells_by_tail[fc_tail]),
                                effect_size_type = "Minimum effect size",
                                tail_type = fc_tail,
                                FDP_eval = sapply(cutoff_grid, function(cutoff_instance){
                                  perturbplan:::compute_FDP_plan(mean_list = max_mean_sd_list$means,
                                                                 sd_list = max_mean_sd_list$sds,
                                                                 side = "left", 
                                                                 prop_non_null = positive_proportion,
                                                                 cutoff = cutoff_instance)}),
                                cutoff_eval = cutoff_grid
                              ),
                              # Median estimate for this tail
                              data.frame(
                                reads_per_cell = unname(reads_per_cell),
                                num_total_cells = unname(num_total_cells),
                                num_trt_cells = unname(num_trt_cells_by_tail[fc_tail]), 
                                effect_size_type = "Median effect size",
                                tail_type = fc_tail,
                                FDP_eval = sapply(cutoff_grid, function(cutoff_instance){
                                  perturbplan:::compute_FDP_plan(mean_list = median_mean_sd_list$means,
                                                                 sd_list = median_mean_sd_list$sds,
                                                                 side = "left", 
                                                                 prop_non_null = positive_proportion,
                                                                 cutoff = cutoff_instance)}),
                                cutoff_eval = cutoff_grid
                              ))
  }
  
}

##################### Plotting the FDP curve ###################################
# select two setups
FDP_versus_cutoff_plot <- analytic_results |> 
  dplyr::group_by(tail_type, cutoff_eval) |> 
  dplyr::filter(num_total_cells == 37888, reads_per_cell %in% c(12287, 22117)) |> 
  dplyr::mutate(reads_per_cell = sprintf("Reads / cell = %d", reads_per_cell),
                num_total_cells = sprintf("Cells = %d", num_total_cells),
                cell_reads = sprintf("%s, %s", reads_per_cell, num_total_cells)) |>
  dplyr::mutate( 
    tail_type = ifelse(tail_type == "light_tail", "Lower effect variability", "Higher effect variability"),
    tail_type = factor(tail_type, levels = c("Lower effect variability", "Higher effect variability"))) |>
  ggplot(aes(x = cutoff_eval, y = FDP_eval, color = effect_size_type)) + 
  geom_line(size = 0.5) + 
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_wrap(cell_reads ~ tail_type, scales = "free_x",
             labeller = labeller(tail_type = function(x) str_replace_all(str_to_title(x), "_", " "))) +  # Same as power plots
  scale_color_manual(
    values = c("Median effect size" = "#fee090", "Minimum effect size" = "#d73027")
  ) +
  labs(x = "Significance cutoff", y = "False discovery proportion (FDP) estimate") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14)
  )
  
# save the plot 
path_to_save <- "reproduction-code-prep/final_plots/figures/"
ggsave(sprintf("%s/FDP_curve_estimation.pdf", path_to_save), FDP_versus_cutoff_plot, width = 7.5, height = 7.5)

