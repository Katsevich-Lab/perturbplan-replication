# This is a Rscript on using PerturbPlan to estimate power for fold change distribution analysis

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
experiment <- args[1]

# specify positive proportion and significance level
if(experiment == "Gasperini"){
  positive_proportion <- 0.01
  alpha <- 0.01
  window_def <- "gene_analysis"
}else{
  positive_proportion <- 0.1
  alpha <- 0.1
  window_def <- "enhancer_analysis"
}

# specify guides per target
if(experiment == "Ray"){
  guides_per_target <- 15
}else{
  guides_per_target <- 2
}

# source the necessary files
source("~/.Rprofile")
source("reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R")
source("reproduction-code-prep/realdata-validation-pipeline/run_fc_distr_construction.R")

# specify intermediate files
intermediate_files_folder <- "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/"

# specify the fc_tail_list (replacing fc_spread_list)
fc_tail_list <- c("light_tail", "medium_tail", "heavy_tail")

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
    
    # extract downsampled guides per target
    if(experiment == "Gasperini"){
      downsampled_guides_per_target <- 2
    }else{
      downsampled_guides_per_target <- unique(plan_parameter_grid$downsampled_guides_per_target)
    }
    
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
    
    ############################### 3.3 Perform PerturbPlan for both approaches #####################
    # Regular approach
    overall_power_max[1, fc_tail] <- perturbplan::compute_power_plan_overall(
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_cntrl_cells,
      library_size = library_size,
      multiple_testing_alpha = alpha,
      prop_non_null = positive_proportion,
      fc_expression_df = fc_expression_df_max$fc_expression_df
    )
    overall_power_median[1, fc_tail] <- perturbplan::compute_power_plan_overall(
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_cntrl_cells,
      library_size = library_size,
      multiple_testing_alpha = alpha,
      prop_non_null = positive_proportion,
      fc_expression_df = fc_expression_df_median$fc_expression_df
    )
  }
  
  # Append both results with tail-specific num_trt_cells in long format
  for (tail in fc_tail_list) {
    analytic_results <- rbind(analytic_results,
                              # Lower bound for this tail
                              data.frame(
                                reads_per_cell = reads_per_cell,
                                num_total_cells = num_total_cells, 
                                num_trt_cells = num_trt_cells_by_tail[tail],
                                power_type = "lower_bound",
                                tail_type = tail,
                                power = overall_power_max[1, tail]
                              ),
                              # Median estimate for this tail
                              data.frame(
                                reads_per_cell = reads_per_cell,
                                num_total_cells = num_total_cells,
                                num_trt_cells = num_trt_cells_by_tail[tail], 
                                power_type = "median_estimate",
                                tail_type = tail,
                                power = overall_power_median[1, tail]
                              ))
  }
  
}

# save the results
path_to_save <- "reproduction-code-prep/realdata-validation-pipeline/PerturbPlan/results"
if(!dir.exists(path_to_save)){
  dir.create(path_to_save, recursive = TRUE)
}
saveRDS(analytic_results, sprintf("%s/%s_%s_perturbplan_overall_power.rds", path_to_save, experiment, window_def))
