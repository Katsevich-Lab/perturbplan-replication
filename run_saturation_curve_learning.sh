#!/bin/bash

#$ -N run_saturation_curve       # Specify job name
#$ -l m_mem_free=80G             # to request 64 GB of memory per core.
#$ -cwd                          # to run the job in the current working directory.

Rscript reproduction-code-prep/saturation-curve-fitting/Gasperini-10x-learning.R
Rscript reproduction-code-prep/saturation-curve-fitting/Gasperini-10x-plotting.R
