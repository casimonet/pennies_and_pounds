#!/bin/bash
#$ -N PP_simulations           
#$ -o .//logs
#$ -e ./logs      
#$ -l h_rt=23:59:00 
#$ -l h_vmem=10G


# Initialise the environment modules
. /etc/profile.d/modules.sh


# for each of these samples, launch the R script
module load R
Rscript --vanilla simulations.R --simulations 20 --seed $SGE_TASK_ID --outdir ./output