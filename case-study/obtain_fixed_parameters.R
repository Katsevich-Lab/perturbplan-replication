# specify fixed parameters
fixed_parameters <- list(
  # fdr level
  fdr_target = 0.1,
  
  # precision level for power search
  precision = 0.005,
  
  # specify mapping efficiency
  mapping_efficiency = 0.704,
  mapping_efficiency_tap = 0.371,
  
  # specify gRNA variability
  gRNA_variability = 0.13,
  
  # default TPM threshold
  tpm_threshold_default = 10,
  
  # default non-null proportion
  prop_non_null = 0.005,
  
  # default side
  side = "left",
  
  # default control group type
  control_group = "complement",
  
  # default parameter grid size
  grid_size = 100, 
  
  # default cost per captured cell
  cost_per_captured_cell = 0.086,
  
  # default cost per million reads
  cost_per_million_reads =  0.374,
  
  # specify power target
  power_target = 0.8,
  
  # specify minimum fold change of interest
  minimum_fold_change = 0.75,
  
  # specify the parameter grid for restrospective analysis
  minimum_fold_change_retrospective_list = seq(0.3, 0.95, length.out = 100),
  tpm_threshold_retrospective_list = 10^{seq(0, 3, length.out = 100)},
  
  # specify the parameter grid for prospective analyses
  tpm_threshold_prospective_list = c(3, 10, 50, 200, 500),
  minimum_fold_change_prospective_list = seq(0.6, 0.8, length.out = 5)
)

# save the fixed parameters
saveRDS(fixed_parameters, "reproduction-code-prep/case-study/fixed_parameters.rds")