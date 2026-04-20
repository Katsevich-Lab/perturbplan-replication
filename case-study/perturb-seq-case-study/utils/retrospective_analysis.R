set.seed(1)

################### Power over minimum detectable effect size ##################
# specify the path to the specific dataframe
path_to_power_over_fc_df <- paste0(dir_to_results, "power_over_fc_df.rds")

# compute the power if power_over_fc_df is not there
if(!file.exists(path_to_power_over_fc_df)){
  # compute the power using Perturb-seq
  power_over_fc_df <- compute_power_plan(
    TPM_threshold = 0,
    minimum_fold_change = minimum_fold_change_retrospective_list,
    cells_per_target = cells_per_target,
    sequenced_reads_per_cell = sequenced_reads_per_cell,
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
    min_power_threshold = max(power_target - 0.1, 0.05),
    max_power_threshold = min(power_target + 0.1, 0.95),
    mapping_efficiency = mapping_efficiency
  ) 
  
  # save the dataframe
  saveRDS(power_over_fc_df, path_to_power_over_fc_df)
}else{
  power_over_fc_df <- readRDS(path_to_power_over_fc_df)
}

###################### Power over TPM analysis threshold #######################
# specify the path to the specific dataframe
path_to_power_over_tpm_df <- paste0(dir_to_results, "power_over_tpm_df.rds")

# compute the power if the dataframe is not there
if(!file.exists(path_to_power_over_tpm_df)){
  # compute the power using Perturb-seq
  power_over_tpm_df <- compute_power_plan(
    TPM_threshold = tpm_threshold_retrospective_list,
    minimum_fold_change = minimum_fold_change,
    cells_per_target = cells_per_target,
    sequenced_reads_per_cell = sequenced_reads_per_cell,
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
    min_power_threshold = max(power_target - 0.1, 0.05),
    max_power_threshold = min(power_target + 0.1, 0.95),
    mapping_efficiency = mapping_efficiency
  ) 
  
  # save the dataframe
  saveRDS(power_over_tpm_df, path_to_power_over_tpm_df)
}else{
  power_over_tpm_df <- readRDS(path_to_power_over_tpm_df)
}

################## tpm_threshold versus minimum effect size ####################
# specify the path to the specific dataframe
path_to_fixed_power_df <- paste0(dir_to_results, "fixing_power_varying_tpm_fc_df.rds")

# compute the power if the dataframe is not available
if(!file.exists(path_to_fixed_power_df)){
  # compute the power using Perturb-seq
  fixed_power_df <- compute_power_plan(
    TPM_threshold = tpm_threshold_retrospective_list,
    minimum_fold_change = minimum_fold_change_retrospective_list,
    cells_per_target = cells_per_target,
    sequenced_reads_per_cell = sequenced_reads_per_cell,
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
    min_power_threshold = max(power_target - 0.1, 0.05),
    max_power_threshold = min(power_target + 0.1, 0.95),
    mapping_efficiency = mapping_efficiency
  ) 

  # combine power dataframe
  varying_fc_tpm_df <- fixed_power_df |>
    dplyr::mutate(minimum_fold_change = round((1 - minimum_fold_change) * 100, 1)) |>
    dplyr::filter((overall_power > power_target - precision) & (overall_power < power_target + precision)) |>
    dplyr::group_by(minimum_fold_change) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup() |>
    dplyr::group_by(TPM_threshold) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup() |>
    as.data.frame()

  # save the dataframe
  saveRDS(varying_fc_tpm_df, path_to_fixed_power_df)
}else{
  varying_fc_tpm_df <- readRDS(path_to_fixed_power_df)
}
