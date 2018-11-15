# pennies_and_pounds
Scripts for Pennies and Pounds paper

# models.R
models functions used in deSolve in the simulations script

# simulations.R
Script for conducting simulations. Designed to conduct simulations in batch using SGE job submission script.
requires 'models.R' to run.
requires a seed to be set, for reproducibility of results at parameters draw step. 
Output:
sims_nb_x_seed_x.RData : R session that contains all deSolve output, to use for dynamics plots
out_data_x_seed_x.txt: table with all public health measures from each simulation, for each model and each type of intervention, parameters, initial states, and indication if numerical issues were detected during the run

