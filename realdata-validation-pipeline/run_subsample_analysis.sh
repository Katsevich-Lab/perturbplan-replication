######################### 0. Set default arguments #############################
max_gb=12          # specify the preferred number of GB for each simulation process
max_hours=2
profile="standard" # HPC profile
B_check=1          # Number of replicates to use for benchmarking (default 3)
benchmark_memory=0 # Logic value; if 0, ignore the max_gb computation

# read command line arguments
while [ $# -gt 0 ]; do
    if [[ $1 == "--"* ]]; then
        v="${1/--/}"
        declare "$v"="$2"
        shift
    fi
    shift
done

# check that required arguments are given
if [ -z "$B" ]; then
  echo "Error: The simulation name needs to be given via the\
 command-line argument B, number of subsample."
  exit
fi
if [ -z "$sim_spec_dir" ]; then
  echo "Error: The simulation specifier directory needs to be given via the\
 command-line argument sim_spec_dir."
  exit
fi
if [ -z "$experiment" ]; then
  echo "Error: The simulation specifier directory needs to be given via the\
 command-line argument experiment."
  exit
fi
if [ -z "$result_dir" ]; then
  echo "Error: The simulation specifier directory needs to be given via the\
 command-line argument result_dir."
  exit
fi
if [ -z "$output_filename" ]; then
  echo "Error: The simulation specifier directory needs to be given via the\
 command-line argument output_filename."
  exit
fi

############################# 1. Set up R packages #############################
# extract the R library path from renv
Rscript -e 'write(paste0(.libPaths(), collapse = ":"), file = "temp.txt")'
RENV_R_LIBS_USER=$(cat temp.txt)
rm temp.txt

# export the R library path
export R_LIBS_USER=$RENV_R_LIBS_USER

# write the R library path into nextflow config
echo "env.R_LIBS_USER = \"$RENV_R_LIBS_USER\"" > nextflow.config

############################## 2. Run the simulation ###########################
# create and extract the simulatr specifier
Rscript reproduction-code-prep/realdata-validation-pipeline/sim-spec/real_data.R "${experiment}" "$B" "${grna_threshold}" "${downsampled_guides_per_target}"
sim_spec_fp=$(realpath "reproduction-code-prep/realdata-validation-pipeline/sim-spec/real_data_${experiment}.rds")

# running sceptre on subsampled data
echo "Running the plan for "${experiment}" experiment ..."
nextflow pull katsevich-lab/simulatr-pipeline -r memory-efficient
nextflow run katsevich-lab/simulatr-pipeline -r memory-efficient \
  --simulatr_specifier_fp $sim_spec_fp \
  --result_dir $result_dir \
  --result_file_name $output_filename \
  --B_check $B_check \
  --B $B \
  --max_gb $max_gb \
  --max_hours $max_hours \
  --benchmark_memory $benchmark_memory \
  -profile $profile \
  -with-trace ${result_dir}/real_data_${experiment}_plan_trace.txt
if [ -f "$result_dir/$output_filename" ]; then
  nextflow clean
fi
echo -e "Real data validation with sceptre for ${experiment} complete.\nVery Good!\n....\n....\n"

