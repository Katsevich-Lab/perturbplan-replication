# This is a script trying to analyze the Morris data
library(dplyr)

# source helper script
source("reproduction-code-prep/case-study/helper.R")

# load K562 pilot data
data("K562_Gasperini", package = "perturbplan")

# load the discovery results
discovery_dataset <- "Replogle"
discovery_table <- obtain_discovery_table(dataset_name = discovery_dataset)

# obtain auxiliary information
replogle_info <- readxl::read_excel("reproduction-code-prep/case-study/discovery-table/replogle_2022_p_values.xlsx", sheet = "TabA_K562_day8_summary_stat")

# path to power results
dir_to_results <- sprintf("reproduction-code-prep/case-study/perturb-seq-case-study/results/%s/", discovery_dataset)
if(!dir.exists(dir_to_results)){
  dir.create(dir_to_results, recursive = TRUE)
}

############################## Obtain fixed parameters #########################
fixed_parameters <- readRDS("reproduction-code-prep/case-study/fixed_parameters.rds")
fdr_target <- fixed_parameters$fdr_target
precision <- fixed_parameters$precision
mapping_efficiency <- fixed_parameters$mapping_efficiency
gRNA_variability <- fixed_parameters$gRNA_variability
tpm_threshold_default <- fixed_parameters$tpm_threshold_default
prop_non_null <- fixed_parameters$prop_non_null
control_group <- fixed_parameters$control_group
side <- fixed_parameters$side
grid_size <- fixed_parameters$grid_size
cost_per_captured_cell <- fixed_parameters$cost_per_captured_cell
cost_per_million_reads <- fixed_parameters$cost_per_million_reads
power_target <- fixed_parameters$power_target
minimum_fold_change <- fixed_parameters$minimum_fold_change

# specify the parameter grid for retrospective analyses
tpm_threshold_retrospective_list <- fixed_parameters$tpm_threshold_retrospective_list
minimum_fold_change_retrospective_list <- fixed_parameters$minimum_fold_change_retrospective_list

# specify the parameter grid for prospective analyses
tpm_threshold_prospective_list <- fixed_parameters$tpm_threshold_prospective_list
minimum_fold_change_prospective_list <- fixed_parameters$minimum_fold_change_prospective_list

############################# specify necessary parameters #####################
# MOI parameter
MOI <- 1

# specify number of total cells
num_total_cells <- round(replogle_info |>
                           dplyr::summarise(sum(`number of cells (filtered)`, na.rm = TRUE)) |>
                           dplyr::pull())

# extract library parameters
library_parameters <- K562_Gasperini$library_parameters

# reads per cell
sequenced_reads_per_cell <- round({
  median_UMI_per_cell <- 1e4
  reads_per_cell_list <- 10^{seq(log10(1e4), log10(1e5), length.out = 100)}
  reads_per_cell_list[which.min(abs(
    fit_read_UMI_curve_cpp(reads_per_cell = reads_per_cell_list * mapping_efficiency, 
                           rSAC_fn_wrapper = library_parameters) - median_UMI_per_cell
  ))]
} / mapping_efficiency)

# extract MOI, gRNAs per target 
gRNAs_per_target <- 1

# specify number of target
num_targets <- 9867

# specify the number of non-targeting gRNAs
non_targeting_gRNAs <- 585

# extract cells per target
cells_per_target <- round(num_total_cells * MOI * gRNAs_per_target / (num_targets + non_targeting_gRNAs))

########################### Extract pilot data #################################
# extract discovery pairs
discovery_pairs <- discovery_table |> 
  dplyr::select(target_id, response_id) |>
  dplyr::rename(grna_target = target_id) |>
  dplyr::distinct()

# extract baseline expression statistics
baseline_expression_stats <- perturbplan:::extract_expression_info(biological_system = "K562",
                                                                   B = 2e3, 
                                                                   gene_list = discovery_pairs$response_id,
                                                                   TPM_threshold = 0)$expression_df

###################### Perform the retrospective analyses ######################
# specify analyses parameters
control_group <- "nt_cells"

# retrospective analysis
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/retrospective_analysis.R")

# do restrospective plotting
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/plot_retrospective_analysis.R")

######################### Perform cost-power analysis ##########################
# prospective analysis
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/prospective_analysis.R")

# do prospective analysis plot
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/plot_prospective_analysis.R")
