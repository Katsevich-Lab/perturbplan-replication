################### Power over minimum detectable effect size ##################
# specify the path to the specific dataframe
path_to_power_over_fc_df <- paste0(dir_to_results, sprintf("power_over_fc_df_at_min_eff_%d_percent.rds", round(primer_efficiency_threshold * 100)))

# compute the power if power_over_fc_df is not there
if(!file.exists(path_to_power_over_fc_df)){
  
  # obtain the filtered genes
  filtered_genes <- unique(
    baseline_expression_stats_tap |> 
      dplyr::filter(UMI_per_cell >= UMI_threshold_default) |> 
      dplyr::select(response_id) |>
      dplyr::pull()
  )
  
  # compute the power using Perturb-seq
  Perturb_over_fc_df <- compute_power_plan(
    TPM_threshold = 0,
    minimum_fold_change = minimum_fold_change_retrospective_list,
    cells_per_target = cells_per_target,
    sequenced_reads_per_cell = reads_per_cell,
    MOI = MOI,
    num_targets = num_targets,
    non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target,
    gRNA_variability = gRNA_variability,
    control_group = control_group,
    side = side,
    multiple_testing_alpha = fdr_target,
    prop_non_null = prop_non_null,
    baseline_expression_stats = baseline_expression_stats_perturb |> 
      dplyr::filter(!is.na(relative_expression) & response_id %in% filtered_genes),
    library_parameters = library_parameters_perturb,
    min_power_threshold = max(power_target - 0.1, 0.05),
    max_power_threshold = min(power_target + 0.1, 0.95),
    mapping_efficiency = mapping_efficiency_perturb
  )
  
  TAP_over_fc_df <- compute_power_plan(
    TPM_threshold = 0,
    minimum_fold_change = minimum_fold_change_retrospective_list,
    cells_per_target = cells_per_target,
    sequenced_reads_per_cell = reads_per_cell,
    MOI = MOI,
    num_targets = num_targets,
    non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target,
    gRNA_variability = gRNA_variability,
    control_group = control_group,
    side = side,
    multiple_testing_alpha = fdr_target,
    prop_non_null = prop_non_null,
    baseline_expression_stats = baseline_expression_stats_tap |> 
      dplyr::filter(!is.na(relative_expression) & response_id %in% filtered_genes),
    library_parameters = library_parameters_tap,
    min_power_threshold = max(power_target - 0.1, 0.05),
    max_power_threshold = min(power_target + 0.1, 0.95),
    mapping_efficiency = mapping_efficiency_tap
  )
  
  # combine power dataframe
  power_over_fc_df <- Perturb_over_fc_df |> 
    dplyr::mutate(assay = "Perturb-seq") |>
    dplyr::bind_rows(TAP_over_fc_df |> dplyr::mutate(assay = "TAP-seq")) |>
    as.data.frame() 
  
  # save the dataframe
  saveRDS(power_over_fc_df, path_to_power_over_fc_df)
}else{
  power_over_fc_df <- readRDS(path_to_power_over_fc_df)
}

###################### Power over TPM analysis threshold #######################
# specify the path to the specific dataframe
path_to_power_over_tpm_df <- paste0(dir_to_results, sprintf("power_over_tpm_df_at_min_eff_%d_percent.rds", round(primer_efficiency_threshold * 100)))

# compute the power if the dataframe is not there
if(!file.exists(path_to_power_over_tpm_df)){
  
  # create null df
  Perturb_over_tpm_df <- data.frame(NULL)
  TAP_over_tpm_df <- data.frame(NULL)
  
  # Loop over TPM threshold list
  for (UMI_threshold in UMI_threshold_retrospective_list) {
    
    # Filter out genes above tpm_threshold
    filtered_genes <- baseline_expression_stats_tap |> 
      dplyr::filter(UMI_per_cell >= UMI_threshold) |> 
      dplyr::select(response_id) |> 
      dplyr::pull()
    
    # compute the power using Perturb-seq
    Perturb_over_tpm_df <- compute_power_plan(
      TPM_threshold = 0,
      minimum_fold_change = minimum_fold_change,
      cells_per_target = cells_per_target,
      sequenced_reads_per_cell = reads_per_cell,
      MOI = MOI,
      num_targets = num_targets,
      non_targeting_gRNAs = non_targeting_gRNAs,
      gRNAs_per_target = gRNAs_per_target,
      gRNA_variability = gRNA_variability,
      control_group = control_group,
      side = side,
      multiple_testing_alpha = fdr_target,
      prop_non_null = prop_non_null,
      baseline_expression_stats = baseline_expression_stats_perturb |> 
        dplyr::filter(response_id %in% filtered_genes),
      library_parameters = library_parameters_perturb,
      min_power_threshold = max(power_target - 0.1, 0.05),
      max_power_threshold = min(power_target + 0.1, 0.95),
      mapping_efficiency = mapping_efficiency_perturb
    ) |> 
      dplyr::mutate(UMI_threshold = UMI_threshold) |>
      dplyr::select(-TPM_threshold) |>
      dplyr::bind_rows(Perturb_over_tpm_df)
    
    # TAP-seq computation
    TAP_over_tpm_df <- compute_power_plan(
      TPM_threshold = 0,
      minimum_fold_change = minimum_fold_change,
      cells_per_target = cells_per_target,
      sequenced_reads_per_cell = reads_per_cell,
      MOI = MOI,
      num_targets = num_targets,
      non_targeting_gRNAs = non_targeting_gRNAs,
      gRNAs_per_target = gRNAs_per_target,
      gRNA_variability = gRNA_variability,
      control_group = control_group,
      side = side,
      multiple_testing_alpha = fdr_target,
      prop_non_null = prop_non_null,
      baseline_expression_stats = baseline_expression_stats_tap |> 
        dplyr::filter(response_id %in% filtered_genes),
      library_parameters = library_parameters_tap,
      min_power_threshold = max(power_target - 0.1, 0.05),
      max_power_threshold = min(power_target + 0.1, 0.95),
      mapping_efficiency = mapping_efficiency_tap
    ) |>
      dplyr::mutate(UMI_threshold = UMI_threshold) |>
      dplyr::select(-TPM_threshold) |>
      dplyr::bind_rows(TAP_over_tpm_df)
  }
  
  # combine power dataframe
  power_over_tpm_df <- Perturb_over_tpm_df |> 
    dplyr::mutate(assay = "Perturb-seq") |>
    dplyr::bind_rows(TAP_over_tpm_df |> dplyr::mutate(assay = "TAP-seq")) |> 
    as.data.frame() 
  
  # save the dataframe
  saveRDS(power_over_tpm_df, path_to_power_over_tpm_df)
}else{
  power_over_tpm_df <- readRDS(path_to_power_over_tpm_df)
}
