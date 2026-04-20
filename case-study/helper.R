# This is a script crating helper functions for TAP-seq case study
library(perturbplan)
source("~/.Rprofile")

# This is a script loading the discovery table for different datasets
obtain_discovery_table <- function(dataset_name){
  # switch over different dataset 
  switch(dataset_name,
         Ray = {
           readRDS("reproduction-code-prep/case-study/discovery-table/DC_TAP_seq_left_sided_pvalues.rds")
         },
         Gasperini = {
           readRDS(paste0(.get_config_path("LOCAL_PERTURBPLAN_DATA_DIR"), "case-studies/gasperini-discovery-pairs.rds"))
         },
         Morris = {
           readxl::read_excel("reproduction-code-prep/case-study/discovery-table/morris_2023_table_s3.xlsx", sheet = "Table S3F", skip = 2)
         },
         Replogle = {
           readRDS(paste0(.get_config_path("LOCAL_PERTURBPLAN_DATA_DIR"), "case-studies/replogle-discovery-pairs.rds"))
         })
}

# This is a helper script obtaining the baseline expression stats based on primer efficiency filtering
get_baseline_expression <- function(primer_efficiency_threshold, pilot_data){
  
  # load Gasperini and 10x example data
  Gasperini_K562 <- perturbplan:::get_pilot_data_from_package("K562")
  Ray_K562 <- perturbplan:::get_pilot_data_from_package("K562_TAP")
  
  # transform the relative expression to TPM
  representation_target_genes <- mean(Gasperini_K562$baseline_expression_stats |>
                                        dplyr::filter(!is.na(relative_expression)) |>
                                        dplyr::filter(response_id %in% unique(pilot_data$response_id)) |>
                                        dplyr::summarise(
                                          representation_Gasperini = sum(relative_expression)
                                        ) |>
                                        dplyr::pull())

  # extract discovery pairs and compute share-normalized primer efficiency:
  # primer_eff_j = (UMI_TAP_j     / sum_{k in targets} UMI_TAP_k)
  #              / (UMI_Perturb_j / sum_{k in targets} UMI_Perturb_k)
  # Genes are retained when primer_eff_j >= primer_efficiency_threshold
  # (larger = better primer; 0 = keep all genes).
  discovery_pairs <- pilot_data |>
    dplyr::select(grna_target, response_id, relative_expression_at_scale) |>
    dplyr::distinct() |>
    dplyr::rename(TAP_seq = relative_expression_at_scale) |>
    dplyr::left_join(
      Gasperini_K562$baseline_expression_stats |>
        dplyr::filter(!is.na(relative_expression)) |>
        dplyr::select(relative_expression, response_id) |>
        dplyr::rename(Perturb_seq = relative_expression) |>
        dplyr::distinct(),
      by = "response_id"
    ) |>
    dplyr::mutate(
      TAP_UMI_per_cell     = Ray_K562$library_parameters$UMI_per_cell      * TAP_seq,
      Perturb_UMI_per_cell = Gasperini_K562$library_parameters$UMI_per_cell * Perturb_seq
    ) |>
    dplyr::filter(!is.na(TAP_UMI_per_cell), !is.na(Perturb_UMI_per_cell)) |>
    dplyr::mutate(
      tap_share     = TAP_UMI_per_cell     / sum(TAP_UMI_per_cell),
      perturb_share = Perturb_UMI_per_cell / sum(Perturb_UMI_per_cell),
      primer_eff    = tap_share / perturb_share
    ) |>
    dplyr::filter(primer_eff >= primer_efficiency_threshold) |>
    dplyr::select(grna_target, response_id)
  
  # Extract the baseline expression data frame
  baseline_expression_stats_perturb <- Gasperini_K562$baseline_expression_stats|>
    dplyr::filter(response_id %in% discovery_pairs$response_id) |> 
    dplyr::mutate(UMI_per_cell = relative_expression * Gasperini_K562$library_parameters$UMI_per_cell)
  baseline_expression_stats_tap <- Ray_K562$baseline_expression_stats|>
    dplyr::filter(response_id %in% discovery_pairs$response_id) |>
    dplyr::distinct() |>
    dplyr::inner_join(baseline_expression_stats_perturb|>
                        dplyr::select(response_id, relative_expression) |>
                        dplyr::rename(relative_expression_main = relative_expression) |>
                        dplyr::distinct(),
                      by = "response_id") |> 
    dplyr::mutate(UMI_per_cell = relative_expression * Ray_K562$library_parameters$UMI_per_cell)

  # Extract library parameters
  library_parameters_perturb <- Gasperini_K562$library_parameters
  library_parameters_tap <- Ray_K562$library_parameters
  
  # Return the data frame
  return(list(
    baseline_expression_stats_perturb = baseline_expression_stats_perturb |> dplyr::filter(!is.na(relative_expression)),
    baseline_expression_stats_tap = baseline_expression_stats_tap |> dplyr::filter(!is.na(relative_expression)),
    library_parameters_perturb = library_parameters_perturb,
    library_parameters_tap = library_parameters_tap
  ))
}
