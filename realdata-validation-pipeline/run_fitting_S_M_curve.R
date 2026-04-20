library(perturbplan)

######################## preprocess the pilot h5 data ##########################
switch(experiment,
       Gasperini = {
         # extract library parameter from perturbplan
         data("K562_Gasperini")
         fitted_S_M_curve <- K562_Gasperini$library_parameters
         
         # specify targeted ratio
         targeted_genes_reads_ratio <- 1
       },
       Ray = {
         # preprocess the molecule_info.h5 for gene expression
         path_to_QC_data <- obtain_file_path(experiment = "Ray")
         sample_dirs_list <- list.files(path = paste0(path_to_QC_data["raw"], "/cell_ranger/"), pattern = "230327", full.names = TRUE)
         num_files_merged <- 3
         pooled_data <- NULL
         for (sample_dir in sample_dirs_list[1:num_files_merged]) {
           
           # process the raw data
           QC_data <- readRDS(paste0(sample_dir, "/processed/QC_data.rds")) 
           
           # rbind the QCed data
           pooled_data <- rbind(pooled_data, QC_data |> dplyr::mutate(SRR = basename(sample_dir)))
         }
         
         # extract the gene reads 
         QC_data <- QC_data[startsWith(QC_data$response_id, "ENSG"), ]
         
         # obtain reads per cell projected to target genes
         target_panel <- readr::read_csv(sprintf("%s/target_panel.csv", sample_dirs_list[1]), 
                                         skip = 5, col_names = TRUE) |> 
           dplyr::select(gene_id) |> 
           dplyr::distinct() |> 
           dplyr::pull()
         targeted_genes_reads <- QC_data |> dplyr::filter(response_id %in% target_panel)
         targeted_genes_reads_ratio <- sum(targeted_genes_reads$num_reads) / sum(QC_data$num_reads)
         
         ############################## fit the S-M cruve ###############################
         fitted_S_M_curve <- perturbplan:::library_estimation(QC_data = targeted_genes_reads)
       }
)
