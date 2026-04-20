#!/bin/bash

#$ -N run_realdata_validation    # Specify job name
#$ -l m_mem_free=24G             # to request 64 GB of memory per core.
#$ -cwd                          # to run the job in the current working directory.

# specify the parameters
experiment=$1
num_subsample=200
positive_proportion=0.01
num_subsampled_genes=1000
grna_threshold=5
downsampled_guides_per_target=10

# restore the package from renv
Rscript -e 'renv::activate(); renv::restore()'

####################### 1. preprocess the planned experiment ###################
Rscript reproduction-code-prep/realdata-validation-pipeline/preprocessing-files/preprocessing_${experiment}.R
Rscript reproduction-code-prep/realdata-validation-pipeline/run_at_scale_sceptre_analysis.R "${experiment}" "${positive_proportion}" "${num_subsampled_genes}"

# ##################### 2. subsample data and perform sceptre ####################
# create the results directory
output_dir="$LOCAL_PERTURBPLAN_DATA_DIR/realdata-validation/${experiment}" 

# specify the simulatr specifier directory
sim_spec_dir=reproduction-code-prep/realdata-validation-pipeline/sim-spec

# make intermediate-files and subsample subfolder
mkdir -p $SCRATCH_DIR/${experiment}/subsample

# create the output directory
mkdir -p $output_dir

# specify the output file names
output_filename="real_data_results.rds"

# if the output_file is already there, then the simulation will be skipped
if [ ! -f "$output_dir/$output_filename" ]; then
  bash reproduction-code-prep/realdata-validation-pipeline/run_subsample_analysis.sh \
    --B $num_subsample \
    --sim_spec_dir $sim_spec_dir \
    --experiment $experiment \
    --result_dir $output_dir \
    --output_filename $output_filename \
    --grna_threshold $grna_threshold \
    --downsampled_guides_per_target $downsampled_guides_per_target
  wait
else
  echo "$output_filename already exists"
fi

######################## 3. Compute the PerturbPlan power ######################
Rscript reproduction-code-prep/realdata-validation-pipeline/PerturbPlan/Prospective_power_perturbplan.R "${experiment}"

################# 4. Summarize downsampling power and plotting #################
Rscript reproduction-code-prep/realdata-validation-pipeline/PerturbPlan/summarize_validation_results.R "${experiment}"

