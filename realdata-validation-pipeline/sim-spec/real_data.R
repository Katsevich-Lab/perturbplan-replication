# This is a Rscript including the analysis for (pilot, plan) = (Gasperini, Gasperini)
library(simulatr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
experiment <- args[1]
B <- as.numeric(args[2])
grna_threshold <- as.numeric(args[3])
downsampled_guides_per_target <- as.numeric(args[4])

# obtain the intermediate-files path
intermediate_files_path <- sprintf("reproduction-code-prep/realdata-validation-pipeline/intermediate-files/%s", experiment)

# obtain the discovery pairs
discovery_pairs <- readRDS(sprintf("%s/discovery_es_at_scale.rds", intermediate_files_path))

# obtain the summary information for SRR
summary_SRR <- readRDS(sprintf("%s/summary_SRR.rds", intermediate_files_path))

#####################
# 1. Preparing varying parameters
#####################
get_ground_truth <- function(n){
  return(NULL)
}

# specify the parameter grids
switch (experiment,
        Gasperini = {
          # specify the varying parameters
          varying_parameters <- list(
            N = round(seq(3, 12, length.out = 10) * summary_SRR$minimum_num_cells / summary_SRR$num_SRRs) * summary_SRR$num_SRRs,
            R = round(seq(0.2, 1, length.out = 10) * summary_SRR$minimum_num_reads_per_cell)
          )
          
          # specify the fixed parameters
          fixed_parameters <- list(
            N = round(6 * summary_SRR$minimum_num_cells / summary_SRR$num_SRRs) * summary_SRR$num_SRRs,
            R = round(0.6 * summary_SRR$minimum_num_reads_per_cell)
          )
          
          # combine parameters
          parameter_grid <- simulatr::create_param_grid_fractional_factorial(varying_values = varying_parameters,
                                                                             baseline_values = fixed_parameters) |>
            dplyr::mutate(downsampled_guides_per_target = "full") |>
            simulatr::add_ground_truth(get_ground_truth)
        },
        Ray = {
          # specify the varying parameters
          varying_parameters <- list(
            N = round(seq(2, 11, length.out = 10) * summary_SRR$minimum_num_cells / summary_SRR$num_SRRs) * summary_SRR$num_SRRs,
            R = round(10^(seq(-2, 0, length.out = 10)) * summary_SRR$minimum_num_reads_per_cell)
          )
          
          # specify the fixed parameters
          fixed_parameters <- list(
            N = round(10 * summary_SRR$minimum_num_cells / summary_SRR$num_SRRs) * summary_SRR$num_SRRs,
            R = round(1.0 * summary_SRR$minimum_num_reads_per_cell)
          )
          
          # combine parameters
          parameter_grid <- simulatr::create_param_grid_fractional_factorial(varying_values = varying_parameters,
                                                                             baseline_values = fixed_parameters) |>
            dplyr::mutate(downsampled_guides_per_target = downsampled_guides_per_target) |>
            simulatr::add_ground_truth(get_ground_truth)
        }
)

#######################
# 2. Preparing fixed parameters
#######################
# specify fixed parameter
fixed_parameters <- list(
  grna_threshold = grna_threshold,                                             # grna threshold for deciding grna assignment
  discovery_pairs = discovery_pairs,                                           # discovery pairs of interest
  experiment = experiment,                                                     # planning experiment
  B = B,                                                                       # number of subsample for each SRR
  seed = 10                                                                    # set random seed
)

###################
# 3. Generate data
###################
# define data-generating model assuming only one enhancer is involved and balanced experimental setup
generate_subsampled_data <- function(R, N, experiment, discovery_pairs, grna_threshold, downsampled_guides_per_target = 10){
  
  ################### source the helper script #################################
  source(paste0(.get_config_path("LOCAL_CODE_DIR"), "experimental-design/reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R"))
  source(paste0(.get_config_path("LOCAL_CODE_DIR"), "experimental-design/reproduction-code-prep/realdata-validation-pipeline/helper-subsample-analysis.R"))

  ################### obtain subsample with the (R, N) #########################
  subsampled_data <- switch (experiment,
    Gasperini = {
      construct_subsample_Gasperini(plan_experiment = experiment,
                                    R = R, N = N,
                                    gene_of_interest = unique(discovery_pairs$response_id),
                                    grna_of_interest = unique(discovery_pairs$grna_id))
    },
    Ray = {
      sample_summary <- readRDS(paste0(.get_config_path("LOCAL_CODE_DIR"),
                                       "experimental-design/reproduction-code-prep/realdata-validation-pipeline/intermediate-files/Ray/summary_persample.rds"))
      construct_subsample_Ray(plan_experiment = experiment,
                              R = R, N = N,
                              gene_of_interest = unique(discovery_pairs$response_id),
                              grna_of_interest = unique(discovery_pairs$grna_id),
                              grna_threshold = grna_threshold,
                              sample_summary = sample_summary)
    }
  )
  
  ############################ save subsampled data ############################
  # create hashing id
  hashing_id <- stringi::stri_rand_strings(n = 1, length = 20)
  
  # save the response matrix
  saveRDS(subsampled_data$response_matrix,
          paste0(.get_config_path("SCRATCH_DIR"),
                 sprintf("/%s/subsample/R_%d_N_%d_%s_response_matrix.rds", experiment, R, N, hashing_id)))
  
  # save the grna matrix
  saveRDS(subsampled_data$grna_matrix,
          paste0(.get_config_path("SCRATCH_DIR"), 
                 sprintf("/%s/subsample/R_%d_N_%d_%s_grna_matrix.rds", experiment, R, N, hashing_id)))
  
  ######################### construct grna data frame ##########################
  # construct the grna_group_data_frame
  grna_group_data_frame <- switch(experiment,
                                  Gasperini = {
                                    discovery_pairs |>
                                      dplyr::select(grna_id, grna_target) |>
                                      dplyr::distinct() 
                                  },
                                  Ray = {
                                    discovery_pairs |>
                                      dplyr::select(grna_id, grna_target) |>
                                      dplyr::distinct() |> 
                                      dplyr::group_by(grna_target) |> 
                                      dplyr::slice_sample(n = downsampled_guides_per_target) |> 
                                      dplyr::ungroup()
                                  })
  
  # return the output
  return(list(R = R, N = N, hashing_id = hashing_id,
              discovery_pairs = discovery_pairs |> dplyr::select(response_id, grna_target) |> dplyr::distinct(),
              grna_group_data_frame = grna_group_data_frame,
              experiment = experiment, gene_covariate = subsampled_data$gene_covariate))
  
}

# construct generate_data_function
generate_data_function <- simulatr::simulatr_function(
  f = generate_subsampled_data,
  arg_names = formalArgs(generate_subsampled_data),
  loop = TRUE
)

######################
# 4. Method functions
######################
# simulated power function with sceptre with covariate
sceptre_power_spec_f <- simulatr::simulatr_function(f = function(data){
  
  ######################## load necessary matrices #############################
  # load the response matrix using hashing_id
  response_matrix <- readRDS(paste0(.get_config_path("SCRATCH_DIR"),
                                    sprintf("/%s/subsample/R_%d_N_%d_%s_response_matrix.rds",
                                            data$experiment, data$R, data$N, data$hashing_id)))
  # load the grna matrix using hashing_id
  grna_matrix <- readRDS(paste0(.get_config_path("SCRATCH_DIR"),
                                sprintf("/%s/subsample/R_%d_N_%d_%s_grna_matrix.rds",
                                        data$experiment, data$R, data$N, data$hashing_id)))
  
  ####################################### perform sceptre ######################
  # create the sceptre object
  sceptre_object <- sceptre::import_data(response_matrix = response_matrix,
                                         grna_matrix = grna_matrix[data$grna_group_data_frame$grna_id, ],
                                         grna_target_data_frame = data$grna_group_data_frame,
                                         extra_covariates = data$gene_covariate |>
                                           tibble::column_to_rownames(var = "cell_id"),
                                         moi = "high")
  
  # set parameter with QC
  sceptre_object <- sceptre::set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = data$discovery_pairs,
    grna_integration_strategy = "union",
    control_group = "complement",
    resampling_mechanism = "permutations",
    side = "left",
    formula_object = formula(~ log(library_size))
  ) |> sceptre::assign_grnas(method = "thresholding", threshold = 1) 
  
  # run_discovery_analysis for sceptre_object
  sceptre_object_w_QC <- sceptre::run_discovery_analysis(sceptre_object |>
                                                           sceptre::run_qc(
                                                             n_nonzero_trt_thresh = 7,
                                                             n_nonzero_cntrl_thresh = 7,
                                                             response_n_umis_range = c(0, 1),
                                                             response_n_nonzero_range = c(0, 1),
                                                             p_mito_threshold = 1
                                                           ))
  sceptre_object_wo_QC <- sceptre::run_discovery_analysis(sceptre_object |>
                                                            sceptre::run_qc(
                                                              n_nonzero_trt_thresh = 0,
                                                              n_nonzero_cntrl_thresh = 0,
                                                              response_n_umis_range = c(0, 1),
                                                              response_n_nonzero_range = c(0, 1),
                                                              p_mito_threshold = 1
                                                            ))
  
  # select the n_nonzero_trt, n_nonzero_ctrl, p-value, response_id and grna_target
  sceptre_discovery_result_w_QC <- sceptre_object_w_QC@discovery_result |>
    dplyr::mutate(pair = sprintf("%s_%s", grna_target, response_id)) |>
    dplyr::select(pair, p_value)
  sceptre_discovery_result_wo_QC <- sceptre_object_wo_QC@discovery_result |>
    dplyr::mutate(pair = sprintf("%s_%s", grna_target, response_id)) |>
    dplyr::select(pair, p_value)
  
  # create a vector with name being the response_id and value being the p_value
  pvalue_vec_w_QC <- setNames(sceptre_discovery_result_w_QC$p_value, sceptre_discovery_result_w_QC$pair)
  pvalue_vec_wo_QC <- setNames(sceptre_discovery_result_wo_QC$p_value, sceptre_discovery_result_wo_QC$pair)
  
  # return the power function
  return(list(
    p_value_w_QC = pvalue_vec_w_QC,
    p_value_wo_QC = pvalue_vec_wo_QC
  ))
  
}, arg_names = character(0), loop = T)

# assemble method function
run_method_functions <- list(sceptre = sceptre_power_spec_f)

##########################
# 5. Check simulation
##########################
simulatr_spec <- simulatr_specifier(
  parameter_grid,
  fixed_parameters,
  generate_data_function,
  run_method_functions
)

##########################
# 6. Save files
##########################
# save the parameter grid
saveRDS(parameter_grid |> dplyr::mutate(grid_id = 1:dplyr::n()), sprintf("%s/real_data_parameter_grid.rds", intermediate_files_path))

# save the simulatr specifier
saveRDS(simulatr_spec, sprintf("reproduction-code-prep/realdata-validation-pipeline/sim-spec/real_data_%s.rds", experiment))

# # run the simulation
# B_in = 1
# results <- check_simulatr_specifier_object(simulatr_spec, B_in = B_in)


