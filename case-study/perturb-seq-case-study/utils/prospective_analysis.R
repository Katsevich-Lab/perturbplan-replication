set.seed(1)

############################ Varying TPM analysis threshold ####################
path_to_opt_tpm_power_df <- paste0(dir_to_results, "power_opt_tpm_df.rds")
path_to_opt_tpm_cost_df <- paste0(dir_to_results, "cost_opt_tpm_df.rds")

# calculate the power if the power dataframe is not available
if(!file.exists(path_to_opt_tpm_power_df) | !file.exists(path_to_opt_tpm_cost_df)){
  power_opt_tpm_df <- perturbplan:::cost_power_computation(
    minimizing_variable = "TPM_threshold",
    fixed_variable = list(TPM_threshold = tpm_threshold_prospective_list, minimum_fold_change = minimum_fold_change),
    MOI = MOI,
    num_targets = num_targets,
    non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target,
    gRNA_variability = gRNA_variability,
    control_group = control_group,
    side = side,
    multiple_testing_alpha = fdr_target,
    prop_non_null = prop_non_null,
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters,
    grid_size = grid_size,
    mapping_efficiency = mapping_efficiency,
    cost_per_captured_cell = cost_per_captured_cell,
    cost_per_million_reads = cost_per_million_reads
  )
  
  # find the optimal design
  optimal_design_list <- perturbplan:::find_optimal_cost_design(cost_power_df = power_opt_tpm_df,
                                                                minimizing_variable = "TPM_threshold",
                                                                power_target = power_target,
                                                                power_precision = precision,
                                                                MOI = MOI,
                                                                num_targets = num_targets,
                                                                non_targeting_gRNAs = non_targeting_gRNAs,
                                                                gRNAs_per_target = gRNAs_per_target,
                                                                cost_per_captured_cell = cost_per_captured_cell,
                                                                cost_per_million_reads = cost_per_million_reads,
                                                                cost_grid_size = 2 * grid_size)
  
  # merge optimal cost and power dataframe
  cost_power_varying_tpm_df <- optimal_design_list$optimal_cost_power_df |> 
    dplyr::mutate(TPM_threshold = factor(TPM_threshold, levels = as.character(tpm_threshold_prospective_list)))
  
  # rename the reads_per_cell column in cost_grid_opt_fc
  cost_grid_opt_tpm <- optimal_design_list$optimal_cost_grid |> 
    dplyr::mutate(TPM_threshold = factor(TPM_threshold, levels = as.character(tpm_threshold_prospective_list)))
  
  # save the dataframe
  saveRDS(cost_power_varying_tpm_df, path_to_opt_tpm_power_df)
  saveRDS(cost_grid_opt_tpm, path_to_opt_tpm_cost_df)
}else{
  cost_power_varying_tpm_df <- readRDS(path_to_opt_tpm_power_df)
  cost_grid_opt_tpm <- readRDS(path_to_opt_tpm_cost_df)
}

############################ Varying Fold change ###############################
# specify path to directory
path_to_opt_fc_power_df <- paste0(dir_to_results, "power_opt_fc_df.rds")
path_to_opt_fc_cost_df <- paste0(dir_to_results, "cost_opt_fc_df.rds")

# compute the power if the power or cost dataframe is not available
if(!file.exists(path_to_opt_fc_power_df) | !file.exists(path_to_opt_fc_cost_df)){
  
  # compte perturb-seq
  power_opt_fc_df <- perturbplan:::cost_power_computation(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(TPM_threshold = tpm_threshold_default, minimum_fold_change = minimum_fold_change_prospective_list),
    MOI = MOI,
    num_targets = num_targets,
    non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target,
    gRNA_variability = gRNA_variability,
    control_group = control_group,
    side = side,
    multiple_testing_alpha = fdr_target,
    prop_non_null = prop_non_null,
    baseline_expression_stats = baseline_expression_stats |> dplyr::filter(!is.na(relative_expression)),
    library_parameters = library_parameters,
    grid_size = grid_size,
    mapping_efficiency = mapping_efficiency,
    cost_per_captured_cell = cost_per_captured_cell,
    cost_per_million_reads = cost_per_million_reads
  )
  
  # find the optimal design
  optimal_design_list <- perturbplan:::find_optimal_cost_design(
    cost_power_df = power_opt_fc_df,
    minimizing_variable = "minimum_fold_change",
    power_target = power_target,
    power_precision = precision,
    MOI = MOI,
    num_targets = num_targets,
    non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target,
    cost_per_captured_cell = cost_per_captured_cell,
    cost_per_million_reads = cost_per_million_reads,
    cost_grid_size = 2 * grid_size
  )
  
  # merge optimal cost and power dataframe
  cost_power_varying_fc_df <- optimal_design_list$optimal_cost_power_df |> 
    dplyr::mutate(minimum_fold_change = factor((1-minimum_fold_change) * 100, levels = as.character((1 - minimum_fold_change_prospective_list) * 100))) 
  
  # rename the reads_per_cell column in cost_grid_opt_fc
  cost_grid_opt_fc <- optimal_design_list$optimal_cost_grid |> 
    dplyr::mutate(minimum_fold_change = factor((1-minimum_fold_change) * 100, levels = as.character((1 - minimum_fold_change_prospective_list) * 100))) 
  
  # save the dataframe
  saveRDS(cost_power_varying_fc_df, path_to_opt_fc_power_df)
  saveRDS(cost_grid_opt_fc, path_to_opt_fc_cost_df)
}else{
  cost_power_varying_fc_df <- readRDS(path_to_opt_fc_power_df)
  cost_grid_opt_fc <- readRDS(path_to_opt_fc_cost_df)
}