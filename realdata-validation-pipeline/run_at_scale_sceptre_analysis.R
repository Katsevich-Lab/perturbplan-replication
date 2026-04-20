# This is a Rscript doing sceptre analysis on at_scale dataset

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
experiment <- args[1]
positive_proportion <- as.numeric(args[2])
num_subsampled_genes <- as.numeric(args[3])

# source the useful helper script
source("~/.Rprofile")
source("reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R")

# specify the path to save the results
path_to_save <- sprintf("reproduction-code-prep/realdata-validation-pipeline/intermediate-files/%s", experiment)

# skip the analysis if the results are already obtained
if(!file.exists(sprintf("%s/discovery_es_at_scale.rds", path_to_save))){
  
  # create the directory
  if(!dir.exists(path_to_save)) dir.create(path_to_save, recursive = TRUE)
  
  ############################### 1. Obtain necessary information ##############
  # load the discovery pairs
  discovery_pairs <- readRDS(sprintf("reproduction-code-prep/realdata-validation-pipeline/intermediate-files/%s/discovery_pairs.rds", experiment))
  gene_list <- unique(discovery_pairs$response_id)
  
  # obtain the data directory
  path_to_processed_data <- obtain_file_path(experiment = experiment)
  
  # depending on the assay, the access to data will be different
  switch(experiment,
         Gasperini = {
           
           # sample discovery pairs to run sceptre analysis
           num_alternative <- nrow(discovery_pairs |> dplyr::filter(target_type %in% c("known_enhancer", "gene_tss")))
           num_null <- num_alternative * (1 - positive_proportion) / positive_proportion
           withr::with_seed(seed = 1, 
                            {
                              subsampled_gene <- sample(gene_list, pmin(num_subsampled_genes, length(gene_list)))
                              filtered_negative_controls <- discovery_pairs |> 
                                dplyr::filter(target_type == "non-targeting" & response_id %in% subsampled_gene) |>
                                dplyr::sample_n(num_null)
                            })
           discovery_filtered <- rbind(filtered_negative_controls, 
                                       discovery_pairs |> dplyr::filter(target_type %in% c("known_enhancer", "gene_tss"))) 
           
           # import the cell covariates, grna and gene expression odm matrices
           cell_covariates <- readRDS(paste0(path_to_processed_data["processed"], "/intermediate/cell_covariates.rds"))
           grna_assignment_odm <- ondisc::read_odm(odm_fp = paste0(path_to_processed_data["processed"], "/processed/grna_assignment/matrix.odm"), 
                                                   metadata_fp = paste0(path_to_processed_data["processed"], "/processed/grna_assignment/metadata.rds"))
           gene_odm <- ondisc::read_odm(odm_fp = paste0(path_to_processed_data["processed"], "/processed/gene/matrix.odm"), 
                                        metadata_fp = paste0(path_to_processed_data["processed"], "/processed/gene/metadata.rds"))
           multimodal_odm <- ondisc::multimodal_ondisc_matrix(covariate_ondisc_matrix_list = list(grna = grna_assignment_odm, gene = gene_odm))
           
           # remove cells with 0 grna expression or gene expression
           ok_cells_grna <- (grna_assignment_odm |> ondisc::get_cell_covariates() |> dplyr::pull(n_nonzero)) != 0L
           ok_cells_gene <- (gene_odm |> ondisc::get_cell_covariates() |> dplyr::pull(n_umis)) != 0L
           multimodal_odm_sub <- multimodal_odm[, ok_cells_grna & ok_cells_gene]
           
           # obtain the grna and gene matrices
           grna_matrix <- multimodal_odm_sub@modalities$grna[[seq(1, nrow(multimodal_odm_sub@modalities$grna)), ]]
           rownames(grna_matrix) <- multimodal_odm_sub@modalities$grna |> ondisc::get_feature_ids()
           response_matrix <- multimodal_odm_sub@modalities$gene[[unique(discovery_filtered$response_id), ]]
           rownames(response_matrix) <- unique(discovery_filtered$response_id)
           
           # set the threshold
           grna_threshold <- 1
         },
         Ray = {
           
           # set the discovery_filtered
           discovery_filtered <- discovery_pairs
           
           # load the sceptre object
           sceptre_object <- readRDS(paste0(path_to_processed_data["raw"], "/sceptre_output/final_sceptre_object.rds"))
           grna_threshold <- sceptre_object@grna_assignment_hyperparameters$threshold
           
           # obtain the response_matrix and grna_matrix
           response_matrix <- sceptre_object@response_matrix[[1]]
           grna_matrix <- sceptre_object@grna_matrix[[1]]
         })
  
  ########################### 2. Perform sceptre analysis ######################
  # import data
  sceptre_object <- sceptre::import_data(response_matrix = response_matrix,
                                         grna_matrix = grna_matrix,
                                         grna_target_data_frame = discovery_filtered |> 
                                           dplyr::select(grna_id, grna_target) |> 
                                           dplyr::distinct(),
                                         moi = "high")
  
  # set analysis parameters and run qc
  sceptre_object <- sceptre::set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_filtered |> dplyr::select(grna_target, response_id) |> dplyr::distinct(),
    grna_integration_strategy = "singleton",
    control_group = "complement",
    resampling_mechanism = "permutations", 
    side = "left",
    formula_object = formula(~ log(response_n_umis))) |>
    sceptre::assign_grnas(method = "thresholding", threshold = grna_threshold) |>
    sceptre::run_qc(
      n_nonzero_trt_thresh = 7,
      n_nonzero_cntrl_thresh = 7,
      response_n_umis_range = c(0, 1),
      response_n_nonzero_range = c(0, 1),
      p_mito_threshold = 1,
    )
  
  # run_discovery_analysis for sceptre_object
  sceptre_object <- sceptre::run_discovery_analysis(sceptre_object)
  
  # set analysis parameters and run qc
  sceptre_object_positive_cntrl <- sceptre::set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_filtered |> 
      dplyr::filter(target_type %in% c("known_enhancer", "gene_tss")) |>
      dplyr::select(grna_target, response_id) |> dplyr::distinct(),
    grna_integration_strategy = "union",
    control_group = "complement",
    resampling_mechanism = "permutations", 
    side = "left",
    formula_object = formula(~ log(response_n_umis))) |>
    sceptre::assign_grnas(method = "thresholding", threshold = grna_threshold) |>
    sceptre::run_qc(
      n_nonzero_trt_thresh = 7,
      n_nonzero_cntrl_thresh = 7,
      response_n_umis_range = c(0, 1),
      response_n_nonzero_range = c(0, 1),
      p_mito_threshold = 1,
    )
  
  # run_discovery_analysis for sceptre_object
  sceptre_object_positive_cntrl <- sceptre::run_discovery_analysis(sceptre_object_positive_cntrl)
  
  ############################ 4. post-process the obtained results ############
  switch(experiment, 
         Gasperini = {
           # select the variables
           sceptre_discovery_result <- sceptre_object@discovery_result |>
             dplyr::filter(!is.na(fold_change)) |>
             dplyr::select(grna_id, response_id, fold_change, se_fold_change) |> 
             dplyr::rename(fold_change_at_scale = fold_change, se_fold_change_at_scale = se_fold_change) |>
             dplyr::distinct() 
           sceptre_discovery_positive_cntrl <- sceptre_object_positive_cntrl@discovery_result |>
             dplyr::filter(!is.na(fold_change)) |>
             dplyr::select(grna_target, response_id, fold_change, se_fold_change) |> 
             dplyr::rename(fold_change_at_scale = fold_change, se_fold_change_at_scale = se_fold_change) |>
             dplyr::distinct() 
           
           # extract the relative expression
           full_response_matrix <- multimodal_odm_sub@modalities$gene[[1:nrow(multimodal_odm_sub@modalities$gene), ]]
           rownames(full_response_matrix) <- multimodal_odm@modalities$gene |> ondisc::get_feature_ids()
           cell_id <- 1:ncol(multimodal_odm_sub@modalities$gene)
         },
         Ray = {
           
           # rename the variables
           sceptre_discovery_result <- sceptre_object@discovery_result |>
             dplyr::select(grna_id, response_id, fold_change, se_fold_change) |> 
             dplyr::rename(fold_change_at_scale = fold_change, se_fold_change_at_scale = se_fold_change) |>
             dplyr::distinct() 
           sceptre_discovery_positive_cntrl <- sceptre_object_positive_cntrl@discovery_result |>
             dplyr::select(grna_target, response_id, fold_change, se_fold_change) |> 
             dplyr::rename(fold_change_at_scale = fold_change, se_fold_change_at_scale = se_fold_change) |>
             dplyr::distinct() 
           
           # extract the full matrix
           full_response_matrix <- response_matrix
           cell_id <- 1:ncol(response_matrix)
         })
  chunk_list <- split(cell_id, cut(seq_along(cell_id), breaks = 10, labels = FALSE))
  sum_expression <- setNames(numeric(nrow(full_response_matrix)), rownames(full_response_matrix))
  for (chunk in 1:length(chunk_list)) {
    sum_expression <- apply(full_response_matrix[, chunk_list[[chunk]]], 1, sum) + sum_expression
  }
  relative_expression <- sum_expression / sum(sum_expression)
  
  # obtain the summary gene information
  gene_summary_info <- data.frame(
    response_id = rownames(response_matrix)
  ) |> dplyr::mutate(
    relative_expression_at_scale = relative_expression[response_id],
    expression_size_at_scale = unlist(sapply(rownames(response_matrix), 
                                             function(response_id) sceptre_object@response_precomputations[[response_id]]$theta))[response_id]
  )
  
  # join the dataframe from sceptre results and discovery_filtered
  discovery_filtered <- discovery_filtered |> dplyr::mutate(grna_gene = sprintf("%s_%s", grna_id, response_id))
  joined_fold_change <- sceptre_discovery_result |>
    dplyr::left_join(gene_summary_info, "response_id") |>
    dplyr::mutate(grna_gene = sprintf("%s_%s", grna_id, response_id)) |>
    dplyr::filter(grna_gene %in% discovery_filtered$grna_gene) |>
    dplyr::left_join(discovery_filtered |> dplyr::select(grna_gene, grna_target, target_type, 
                                                         num_oracle_cells, num_total_plan_cells), by = "grna_gene")
  
  ########################### 5. Save the results ##############################
  # save the joined dataframe
  saveRDS(joined_fold_change, sprintf("%s/discovery_es_at_scale.rds", path_to_save))
  saveRDS(sceptre_discovery_positive_cntrl, sprintf("%s/discovery_es_at_scale_positive.rds", path_to_save))
}
