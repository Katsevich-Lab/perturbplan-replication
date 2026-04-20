# This is a script trying to analyze the Morris data
library(dplyr)

# source helper script
source("reproduction-code-prep/case-study/helper.R")

# load K562 pilot data
data("K562_Gasperini", package = "perturbplan")

# load the discovery results
discovery_dataset <- "Morris"
discovery_table <- obtain_discovery_table(dataset_name = discovery_dataset)

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
MOI <- 13

# specify number of total cells
num_total_cells <- 46583

# reads per cell
sequenced_reads_per_cell <- round(60000 / mapping_efficiency)

# extract gRNAs per target 
grna_df <- readxl::read_excel("reproduction-code-prep/case-study/discovery-table/morris_2023_table_s3.xlsx", 
                              sheet = "Table S3A", skip = 2)
gRNAs_per_target <- round(
  grna_df |>
    dplyr::group_by(Target) |>
    dplyr::summarise(grna_per_target = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::summarise(mean(grna_per_target)) |>
    dplyr::pull()
)

# specify non-targeting gRNAs
non_targeting_gRNAs <- nrow(grna_df |> dplyr::filter(Target == "nt"))

# specify number of target
num_targets <- length(unique(grna_df |> 
                               dplyr::filter(Target != "nt") |>
                               dplyr::select(Target) |> 
                               dplyr::pull()))

# extract cells per target
cells_per_target <- round(num_total_cells * MOI * gRNAs_per_target / nrow(grna_df))

########################### Extract pilot data #################################
# extract discovery pairs
discovery_pairs <- discovery_table |> 
  dplyr::select(gRNAs, `Ensembl ID`) |>
  dplyr::rename(grna_target = gRNAs, response_id = `Ensembl ID`)
# obtain baseline expression statistics
baseline_expression_stats <- perturbplan:::extract_expression_info(biological_system = "K562",
                                                                   B = 2e3, 
                                                                   gene_list = discovery_pairs$response_id,
                                                                   TPM_threshold = 0)$expression_df
# obtain library parameters
library_parameters <- K562_Gasperini$library_parameters

###################### Perform the retrospective analyses ######################
# retrospective analysis
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/retrospective_analysis.R")

# do retrospective plotting
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/plot_retrospective_analysis.R")

######################### Perform cost-power analysis ##########################
# prospective analysis
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/prospective_analysis.R")

# do prospective analysis plot
source("reproduction-code-prep/case-study/perturb-seq-case-study/utils/plot_prospective_analysis.R")
