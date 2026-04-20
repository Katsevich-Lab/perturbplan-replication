set.seed(1)

#################### Varying Cost ratio and efficiency bound ###################
# specify path to directory
path_to_opt_power_eff_df <- paste0(dir_to_results, sprintf("power_opt_eff_df_at_min_eff_%d_percent.rds", round(primer_efficiency_threshold * 100)))
path_to_param_df <- paste0(dir_to_results, "cost_opt_eff_param_grid.rds")

# construct necessary information
optimal_power_df_eff <- data.frame(NULL)
baseline_ratio <- cost_per_captured_cell / cost_per_million_reads

# specify the varying parameters; primer_efficiency_bound is a *lower bound*
# on the share-normalized TAP/Perturb ratio (0 = keep all genes).
varying_parameters <- list(
  ratio = baseline_ratio * c(1 / 25, 1 / 5, 1, 5, 25),
  primer_efficiency_bound = c(0, 0.1, 0.2, 0.3, 0.4)
)

# specify the fixed parameters (headline = 40% primer efficiency)
fixed_parameters <- list(ratio = baseline_ratio, primer_efficiency_bound = 0.4)

# construct the parameter grid
parameter_grid <- simulatr::create_param_grid_fractional_factorial(varying_values = varying_parameters,
                                                                   baseline_values = fixed_parameters)

# save the parameter grid 
saveRDS(parameter_grid, path_to_param_df)

# compute the power if the power or cost dataframe is not available
if(!file.exists(path_to_opt_power_eff_df)){
  
  # Loop over cost ratio
  for(setup_id in 1:nrow(parameter_grid)){
    
    # extract ratio 
    cost_ratio <- parameter_grid[setup_id, ]$ratio
    primer_efficiency_bound <- parameter_grid[setup_id, ]$primer_efficiency_bound
    
    # Consider varying the ratio of cost_per_cell / cost_per_million reads
    extracted_pilot_data <- get_baseline_expression(primer_efficiency_threshold = primer_efficiency_bound, 
                                                    pilot_data = pilot_data)
    
    # compute perturb-seq
    power_opt_eff_df_perturb <- perturbplan:::cost_power_computation(
      minimizing_variable = "cost",
      fixed_variable = list(TPM_threshold = 0, minimum_fold_change = minimum_fold_change),
      MOI = MOI,
      num_targets = num_targets,
      non_targeting_gRNAs = non_targeting_gRNAs,
      gRNAs_per_target = gRNAs_per_target,
      gRNA_variability = gRNA_variability,
      control_group = control_group,
      side = side,
      multiple_testing_alpha = fdr_target,
      prop_non_null = prop_non_null,
      baseline_expression_stats = extracted_pilot_data$baseline_expression_stats_perturb,
      library_parameters = extracted_pilot_data$library_parameters_perturb,
      grid_size = grid_size,
      mapping_efficiency = mapping_efficiency_perturb,
      cost_per_captured_cell = cost_per_million_reads * cost_ratio,
      cost_per_million_reads = cost_per_million_reads,
      power_range = 1.6
    )
    
    power_opt_eff_df_tap <- perturbplan:::cost_power_computation(
      minimizing_variable = "cost",
      fixed_variable = list(TPM_threshold = 0, minimum_fold_change = minimum_fold_change),
      MOI = MOI,
      num_targets = num_targets,
      non_targeting_gRNAs = non_targeting_gRNAs,
      gRNAs_per_target = gRNAs_per_target,
      gRNA_variability = gRNA_variability,
      control_group = control_group,
      side = side,
      multiple_testing_alpha = fdr_target,
      prop_non_null = prop_non_null,
      baseline_expression_stats = extracted_pilot_data$baseline_expression_stats_tap,
      library_parameters = extracted_pilot_data$library_parameters_tap,
      grid_size = grid_size,
      mapping_efficiency = mapping_efficiency_tap,
      cost_per_captured_cell = cost_per_million_reads * cost_ratio,
      cost_per_million_reads = cost_per_million_reads,
      power_range = 1.6
    )
    
    # merge optimal cost and power dataframe
    cost_power_varying_eff_df <- power_opt_eff_df_perturb |> 
      dplyr::mutate(assay = "Perturb-seq") |>
      dplyr::bind_rows(power_opt_eff_df_tap |> dplyr::mutate(assay = "TAP-seq")) |>
      dplyr::mutate(minimum_fold_change = factor((1 - minimum_fold_change) * 100)) |> 
      dplyr::filter(overall_power > (power_target - precision) & overall_power < (power_target + precision))
    
    # concatenate the dataframe
    optimal_power_df_eff <- optimal_power_df_eff |> 
      dplyr::bind_rows(cost_power_varying_eff_df |> dplyr::mutate(cost_ratio = cost_ratio, primer_efficiency_bound = primer_efficiency_bound))
  }
  
  # save the RDS file
  saveRDS(optimal_power_df_eff, path_to_opt_power_eff_df)
  
}else{
  optimal_power_df_eff <- readRDS(path_to_opt_power_eff_df)
}
