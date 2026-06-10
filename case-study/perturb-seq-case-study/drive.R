# This is a script running the design comparison
set.seed(1)

# Run Morris analysis
source("case-study/perturb-seq-case-study/design_comparison/run_Morris_analysis.R")

# Run Gasperini analysis
source("case-study/perturb-seq-case-study/design_comparison/run_Gasperini_analysis.R")

# Run Replogle analysis
source("case-study/perturb-seq-case-study/design_comparison/run_Replogle_analysis.R")

# Run 10x analysis
source("case-study/perturb-seq-case-study/design_comparison/run_10x_analysis.R")
