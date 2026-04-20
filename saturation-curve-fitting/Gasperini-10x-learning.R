# load path to raw data directory
library(perturbplan)
library(ggplot2)
library(scales)
source("~/.Rprofile")
source("reproduction-code-prep/saturation-curve-fitting/helper.R")

# load learned pilot data
data("K562_Gasperini")

# specify the number of SRR used for learning read-umi curve
num_SRR_in_use <- 8

# specify downsampling ratio
downsampling_ratio_list <- 10^{seq(-2, 0, length.out = 20)}

# make directory
dir_to_results <- "reproduction-code-prep/saturation-curve-fitting/results"
if(!dir.exists(dir_to_results)){
  dir.create(dir_to_results)
}

# load the read_umi_downsampling if it is computed before
if(file.exists(sprintf("%s/read_umi_downsampling_10x_Gasperini.rds", dir_to_results))){
  read_umi_downsampling <- readRDS(sprintf("%s/read_umi_downsampling_10x_Gasperini.rds", dir_to_results))
  library_parameters <- readRDS(sprintf("%s/library_parameters_on_%d_SRRs.rds", dir_to_results, num_SRR_in_use))
}else{
  
  # preprocess the 10x K562 data
  read_umi_df_10x <- obtain_read_umi_table_10x(downsampling_ratio = downsampling_ratio_list)
  
  # learn the Gasperini K562 data with 8 SRRs
  fit_Gasperini_read <- obtain_read_umi_table_Gasperini(downsampling_ratio = downsampling_ratio_list, num_SRR = num_SRR_in_use)
  
  # preprocess the Gasperini K562 data
  read_umi_df_Gasperini <- fit_Gasperini_read$read_umi_df
  
  # combined two dataframes
  read_umi_downsampling <- read_umi_df_10x |> dplyr::bind_rows(read_umi_df_Gasperini)
  library_parameters <- fit_Gasperini_read$library_parameters
  
  # save the read-UMI dataframe
  saveRDS(read_umi_downsampling, sprintf("%s/read_umi_downsampling_10x_Gasperini.rds", dir_to_results))
  saveRDS(library_parameters, sprintf("%s/library_parameters_on_%d_SRRs.rds", dir_to_results, num_SRR_in_use))
}

###################### construct the S-M curve #################################
# construct read_umi_grid
reads_per_cell_grid <- 10^{seq(2, 4.8, length.out = 500)}
read_umi_learning <- data.frame(
  reads_per_cell = reads_per_cell_grid,
  library_size = perturbplan:::fit_read_UMI_curve_cpp(
    reads_per_cell = reads_per_cell_grid,
    rSAC_fn_wrapper = K562_Gasperini$library_parameters
  ),
  fitted_approach = "1 run"
) |>
  dplyr::bind_rows(
    data.frame(
      reads_per_cell = reads_per_cell_grid,
      library_size = perturbplan:::fit_read_UMI_curve_cpp(
        reads_per_cell = reads_per_cell_grid,
        rSAC_fn_wrapper = library_parameters
      ),
      fitted_approach = ifelse(num_SRR_in_use == 1, "1 run", sprintf("%d runs", num_SRR_in_use))
    )
  )

# save the Read-UMI learning dataframe 
saveRDS(read_umi_learning, sprintf("%s/read_umi_learning_10x_Gasperini.rds", dir_to_results))
