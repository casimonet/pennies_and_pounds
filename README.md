# pennies_and_pounds
Scripts for Pennies and Pounds paper

# simulations_output_samples
Set of examples simulations output (RData files). Those are the ones used to produce example dynamics in figure 2c of the main text. Can be used to run do_dynamic_plots.R to visualise more examples of dynamics.

# Dynamic_plots.pdf
Output of do_dynamic_plots.R script. For this run we selected a set of simulations output that produced various dynamics in order to illustrate the variety of possible dynamics.

# do_dynamic_plots.R
R script to produce plot showing the dynamic of all strains in each compartment and intervention type.

# job_submission.sh
SGE job submission script

# models.R
models functions used in deSolve in the simulations script

# simulations.R
Script for conducting simulations. Designed to conduct simulations in batch using SGE job submission script.
Requires 'models.R' to run.
Requires a seed to be set, for reproducibility of results at parameters draw step. 
Output:
sims_nb_x_seed_x.RData : R session containing all deSolve output, do_dynamic_plots.R can be run on those RData outputs directly to produce the dynamics plots. Four examples of these Rdata files are provided.
out_data_x_seed_x.txt: table with all public health measures from each simulation, for each model and each type of intervention + records of parameters used, initial states (i.e. end of burn-in phase), and indication if numerical issues were detected during the run.


