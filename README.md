# pennies_and_pounds
Scripts for Pennies and Pounds paper

# models.R
models functions used in deSolve in the simulations script

# simulations.R
Script for conducting simulations. Designed to conduct simulations in batch using SGE job submission script.
Requires 'models.R' to run.
Requires a seed to be set, for reproducibility of results at parameters draw step. 
Output:
sims_nb_x_seed_x.RData : R session containing all deSolve output, dynamics plots are produced from these output.
out_data_x_seed_x.txt: table with all public health measures from each simulation, for each model and each type of intervention + records of parameters used, initial states (i.e. end of burn-in phase), and indication if numerical issues were detected during the run.

