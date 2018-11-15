###############################################
#                Pound&Pennies paper          #
#                   Simulations               #
#               Last edit 2018/11/15          #
#                Camille Simonet              #
###############################################


#!/usr/bin/env Rscript

# SET UP----

library("optparse")
library(ggplot2)
library(deSolve)
library(gridExtra)



option_list = list(
  make_option(c("-s", "--simulations"), type="integer", default=NULL,
              help="number of simulations to run", metavar="integer"),
  make_option(c("-d", "--seed"), type="integer", default=NULL,
              help="seed for numbers generation", metavar="integer"),
  make_option(c("-o", "--outdir"), type="character", default="/exports/csce/eddie/biology/groups/mcnally/camille/PoundsandPennies/NEW_STUFF/",
              help="path to output directory", metavar="character")
  
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$simulations)){
  print_help(opt_parser)
  stop("provide a number of simulations to run", call.=FALSE)
}
if (is.null(opt$seed)){
  print_help(opt_parser)
  stop("provide a number for the seed. This might be the array number if this job is launched as an array job.", call.=FALSE)
}



source('models.R')
setwd(paste(opt$outdir))


###### SIMULATIONS SUMMARY ######
# draw a parameter combination
# run a burn in phase
# check presence of resistance at the end of burn in phase, if too low, draw a new parameter combination
# store end of burn in phase state for initial state
# runs the different models using that parameter combination and this initial state
# record variables



# INITIAL CONDITIONS AND PARAMETERS DRAW ----
# draw a parameters combination from some defined ranges to use for model run
# each draw contains all parameters, but not all of them is necessarily used by all models, only the releavant ones for the given model
# so all models are run on same sample of parameters +/- the ones specific to each model
# Constraints on parameters combinations draw filtering:
# not rounded so never equal to zero
# bb < bw, because we set bb = [0-1] * bw
# must lead to a minimum level of resistance (defined by ......) at the end of burn in phase
# then for paper figure plots, we filtered simulations to keep only releavant ones: those where front line drug is used more
# This mean we kept only simulations where the proportion of front line drug usage, measured at the end of burn in phase, is greater than 0.5


nb_sim<- opt$simulations  # number of simulations
my_seed<- opt$seed        # set seed for reproducibility of parameters draw


# parameters storage
pars<- data.frame(bw = numeric(), bb = numeric(), c = numeric(), d = numeric(), g = numeric(),
                  z = numeric(), a = numeric(), e = numeric(),  p = numeric(), mu = numeric(), f = numeric(),
                  my_seed = numeric(), draw = numeric()) # 'draw' used to know how many time parameters were resampled for a given seed (in case a parameter combination had to be discarded), for record and reproducibility of results


initial_states<- vector('list', length = nb_sim)  # initial states storage
outs_burn_in<- vector('list', length = nb_sim)    # to store outputs of burn in phase




burn_in_time <- seq(0, 300, by = 0.05) # burn in phase time and step size


state_burn_in = c(C_00_u = 0, C_00_f = 0, I_00_f = 0.5, I_00_l = 0, # initial states, only front line treated individuals, of all resistance profiles possible, but with a total 5% of resistance
                  C_10_u = 0, C_10_f = 0, I_10_f = 0.5/57, I_10_l = 0,  
                  C_01_u = 0, C_01_f = 0, I_01_f = 0.5/57, I_01_l = 0, 
                  C_11_u = 0, C_11_f = 0, I_11_f = 0.5/57, I_11_l = 0,
                  D = 0, dTrans = 0)



i=1
set.seed(my_seed)


# For adjustment of initial conditions, see later
# DIAGNOSTIC MODEL
initial_states_front_diag<- vector('list', length = nb_sim)
initial_states_last_diag<- vector('list', length = nb_sim)

# NEW DRUG MODEL
initial_states_front_new<- vector('list', length = nb_sim)
initial_states_last_new<- vector('list', length = nb_sim)
initial_states_middle_new<- vector('list', length = nb_sim)



# For model run 
times<- seq(0, 600, 0.05)           # time with 0.05 step sizes
records<- seq(0, tail(times,1), 1)  # to record only 300 points


# To store simulations
runs<- vector('list', length = nb_sim)


# To store values of interest
out_data<- data.frame(
  
  # PUBLIC HEALTH MEASURES
  
  tot_prevalence_no_int = numeric(),
  tot_prevalence_adj_f = numeric(),
  tot_prevalence_adj_l = numeric(),
  tot_prevalence_diag_f = numeric(),
  tot_prevalence_diag_l = numeric(),
  tot_prevalence_new_f = numeric(),
  tot_prevalence_new_l = numeric(),
  tot_prevalence_new_m = numeric(),
  
  nb_death_no_int = numeric(),
  nb_death_adj_f = numeric(),
  nb_death_adj_l = numeric(),
  nb_death_diag_f = numeric(),
  nb_death_diag_l = numeric(),
  nb_death_new_f = numeric(),
  nb_death_new_l = numeric(),
  nb_death_new_m = numeric(),
  
  total_trans_no_int = numeric(),
  total_trans_adj_f = numeric(),
  total_trans_adj_l = numeric(),
  total_trans_diag_f = numeric(),
  total_trans_diag_l = numeric(),
  total_trans_new_f = numeric(),
  total_trans_new_l = numeric(),
  total_trans_new_m = numeric(),
  
  infection_lenght_no_int = numeric(),
  infection_lenght_adj_f = numeric(),
  infection_lenght_adj_l = numeric(),
  infection_lenght_diag_f = numeric(),
  infection_lenght_diag_l = numeric(),
  infection_lenght_new_f = numeric(),
  infection_lenght_new_l = numeric(),
  infection_lenght_new_m = numeric(),
  
  proba_death_no_int = numeric(),
  proba_death_adj_f = numeric(),
  proba_death_adj_l = numeric(),
  proba_death_diag_f = numeric(),
  proba_death_diag_l = numeric(),
  proba_death_new_f = numeric(),
  proba_death_new_l = numeric(),
  proba_death_new_m = numeric(),
  
  multi_no_int = numeric(),
  multi_adj_f  = numeric(),
  multi_adj_l = numeric(),
  multi_diag_f = numeric(),
  multi_diag_l = numeric(),
  multi_new_f = numeric(),
  multi_new_m = numeric(),
  multi_new_l = numeric(),
  
  prop_f_usage = numeric(),
  prop_inapropriate_usage = numeric(),
  
  # PARAMETERS
  
  bw = numeric(),
  bb = numeric(),
  c = numeric(),
  d = numeric(),
  g = numeric(),
  z = numeric(),
  a = numeric(),
  e = numeric(),
  p = numeric(),
  mu = numeric(),
  f = numeric(),
  seed = numeric(),
  draw = numeric(),
  
  # INITIAL CONDITIONS = numeric() END OF BURN IN PHASE,
  
  ini_C_00_u = numeric(),
  ini_C_00_f = numeric(),
  ini_I_00_f = numeric(),
  ini_I_00_l = numeric(),
  ini_C_10_u = numeric(),
  ini_C_10_f = numeric(),
  ini_I_10_f = numeric(),
  ini_I_10_l = numeric(),
  ini_C_01_u = numeric(),
  ini_C_01_f = numeric(),
  ini_I_01_f = numeric(),
  ini_I_01_l = numeric(),
  ini_C_11_u = numeric(),
  ini_C_11_f = numeric(),
  ini_I_11_f = numeric(),
  ini_I_11_l = numeric(),
  ini_D = numeric(),
  ini_dTrans = numeric(),
  
  neg_values_burn_in = numeric(), # to check for integration issues during burn in
  neg_values_run = numeric()      # to check for integration issues during model run
  
)



min_freq_res<- 0.1 # we want to reach minimum 10% of resistance to validate a parameters combination
min_numeric<- 0.1  # minimum absolute value to reach at end of burn in phase to validate a parameter combination (small numbers lead to issues and misleading values of benefits statistics)
index = 0   # to record draw number, index = index + 1 everytime we have to re-draw parameters for a given simulation

while(i < nb_sim + 1){
  
  reached_min_resistance<- FALSE # condition to validate a parameter draw
  reached_min_numeric<- FALSE # condition to validate a parameter draw
  
  while(reached_min_resistance == FALSE | reached_min_numeric == FALSE){
    

    index = index + 1 
    print(paste('burn in', i, index, sep=' : '))
    
      bw<-runif(1,2,5)
      bb<-runif(1,0,1) * bw
      c<-runif(1,0,1)
      d<-runif(1,0,2)
      g<-runif(1,0,2)
      z<-runif(1,0,1)
      a<-rgamma(1, 1, 0.3) 
      e<-runif(1,0,1) * g 
      p<-runif(1,0,1)
      mu<-runif(1,0,0.1)
      f<-runif(1,0,4)
      
      pars_burn<- c(bw=bw, bb=bb, c=c, d=d, g=g, z=z, a=a, e=e, p=p, mu=mu, f=f, my_seed = my_seed, draw = index)
      
      
    # RUN BURN IN PHASE WITH THE SAMPLED PARAMETERS
    outs_burn_in[[i]]<- ode(y = state_burn_in, times = burn_in_time, func = mod1, parms = pars_burn)     # model 1 = no intervention model for this burn in phase
    outs_burn_in[[i]]<- outs_burn_in[[i]][outs_burn_in[[i]][,'time'] %in% seq(0, tail(burn_in_time,1)),] # record only 300 points
    
   
    
    reached_min_resistance<- (sum(tail(outs_burn_in[[i]], 1)[,c('C_10_u', 'C_10_f', 'I_10_f', 'I_10_l', 'C_01_u', 'C_01_f', 'I_01_f', 'I_01_l', 'C_11_u', 'C_11_f', 'I_11_f', 'I_11_l')])/
                             sum(tail(outs_burn_in[[i]], 1)[,c('C_00_u', 'C_00_f', 'I_00_f', 'I_00_l','C_10_u', 'C_10_f', 'I_10_f', 'I_10_l','C_01_u', 'C_01_f', 'I_01_f', 'I_01_l', 'C_11_u', 'C_11_f', 'I_11_f', 'I_11_l')]))>
                             min_freq_res
    
    
    
    reached_min_numeric<- sum(tail(outs_burn_in[[i]], 1)[,c('C_10_u', 'C_10_f', 'I_10_f', 'I_10_l',
                                                            'C_01_u', 'C_01_f', 'I_01_f', 'I_01_l',
                                                            'C_11_u', 'C_11_f', 'I_11_f', 'I_11_l',
                                                            'C_00_u', 'C_00_f', 'I_00_f', 'I_00_l')]) > min_numeric
    
    
  }
  
  
    # ONCE PARAMETERS VALIDATED, RECORD: parameters, initial states, proportion of F usage, amount of bystander selection
    
    pars[i,]<- pars_burn                                  # record parameter combination
    initial_states[[i]]<- tail(outs_burn_in[[i]], 1)[,-1] # record initial states
    
    
    # modify dTrans initial state before the run --> number of infections at the end of burn in phase
    initial_states[[i]]['dTrans'] <- sum(initial_states[[i]][c('I_00_f', 'I_00_l', 'I_10_f', 'I_10_l', 'I_01_f', 'I_01_l', 'I_11_f', 'I_11_l')])



# BEFORE MODEL RUN, WE MUST ADJUST INITIAL CONDITIONS TO MIMIC WHAT THE INTERVENTION IMPLEMENTATION WOULD DO

# DIAGNOSTIC MODEL
    
# in the model, we do not have anymore front line resistants being front line treated in the symptomatic class, i.e. No I_10_f or I_11_f
# But in burn-in phase (before deployment of the diagnostic) we do
# so it IS releavant to keep those classe in the burn in phase
# but at the time of the diagnostic starting, all I_10_f and I_10_f must be converted (i.e. this is the initial implementation of the diagnostic)
# All I_10_f  are transferred to I_10_l, and I_11_f are transferred to I_11_l
# So we modify the initial conditions to have I_10_l = I_10_l + I_10_f , I_11_l = I_11_l + I_11_f 
# then, for code convenience (to avoid dealing with dataframes of different sizes), I keep the I_10_f  and I_11_f  classes but set their rate of change to zero in the model (I_10_f  = 0, I_11_f  = 0), so that they don't matter
# Similarly, for last line diagnostic, I set I_01_f = I_01_f + I_01_l, I_11_f = I_11_f + I_11_l and then I_01_l = 0, I_11_l = 0


# NEW DRUG MODEL
# For the new drug model, there is only one function, describing the dynamic with three drugs F, M, L
# The only difference to model is if we introduce as front, middle or last is the change in initial conditions
# - For introduction as front:
# all resistance to previously F drug become resistant to M
# resistance to F is now zero because it is a new drug
# L remain what there were
    
# - For introduction as last:
# resistance to F remain F
# All resistance to L now are resistant to M
# resistance to L is zero because it is a new drug

# - For introduction as middle:
# it actually does not change anything, those that were resistant to former front line druf remain resistant to F because it is still the same
# and same of those resistant for last line drug

# In addition, we transition from a system of 4*4+1 equation to 5*8+1 equations because of the introduction of the new drug
# So the table of initial conditions must be extended

# NO INTERVENTION MODEL
# No need of adjusting the initial conditions to run this model

# ADJUVANT MODEL
# No need of adjusting the initial conditions to run this model


    
  
  # FOR DIAGNOSTIC MODEL
  # I_10_l = I_10_l + I_10_f , I_11_l = I_11_l + I_11_f , I_10_f  = 0, I_11_f  = 0
  
  
  initial_states_front_diag[[i]] = initial_states[[i]]  # copy the initial states obtained from burn in phase
    
  initial_states_front_diag[[i]][c('I_10_l', 'I_11_l')] =    # modify appropriately
    initial_states_front_diag[[i]][c('I_10_l', 'I_11_l')] + 
    initial_states_front_diag[[i]][c('I_10_f', 'I_11_f')] 
  initial_states_front_diag[[i]][c('I_10_f', 'I_11_f')] = 0
  
  # LAST LINE DIAGNOSTIC
  # I_01_f = I_01_f + I_01_l, I_11_f = I_11_f + I_11_l, I_01_l = 0, I_11_l = 0
  
  initial_states_last_diag[[i]] = initial_states[[i]] 
  
  initial_states_last_diag[[i]][c('I_01_f', 'I_11_f')] = 
    initial_states_last_diag[[i]][c('I_01_f', 'I_11_f')] + 
    initial_states_last_diag[[i]][c('I_01_l', 'I_11_l')] 
  initial_states_last_diag[[i]][c('I_01_l', 'I_11_l')] = 0
  
  
  
  # FOR NEW DRUG MODEL
  # Treatments change with the intervention, those that were treated with front line drug, switch to the new front line drug when this one is added to the market
  # and from now on, follow normal procedure with escalation to middle drug (fomer front line), and then last line (which has not changed)
  # For the resistance profiles, when introducing as front line: 00 -> 000 , 10 -> 010 , 01 -> 001 , 11 -> 011, any other profile is 0 because this front line drug is new so there is no resistance to it so all 100, 110, 101, 111 = 0 
  # basically we add a zero in front of resistance profiles ...
  # and all Xxxx,m = 0, because everyone either switch to the NEW front line drug, or remain with the former last line one
  
  initial_states_front_new[[i]]<-     c(C_000_u = initial_states[[i]][['C_00_u']],
                                        C_000_f = initial_states[[i]][['C_00_f']],
                                        
                                        I_000_f = initial_states[[i]][['I_00_f']],
                                        I_000_m = 0,
                                        I_000_l = initial_states[[i]][['I_00_l']],
                                        
                                        C_100_u = 0,
                                        C_100_f = 0,
                                        
                                        I_100_f = 0,
                                        I_100_m = 0,
                                        I_100_l = 0,
                                        
                                        C_010_u = initial_states[[i]][['C_10_u']],
                                        C_010_f = initial_states[[i]][['C_10_f']],
                                        
                                        I_010_f = initial_states[[i]][['I_10_f']],
                                        I_010_m = 0,
                                        I_010_l = initial_states[[i]][['I_10_l']],
                                        
                                        C_001_u = initial_states[[i]][['C_01_u']],
                                        C_001_f = initial_states[[i]][['C_01_f']],
                                        
                                        I_001_f = initial_states[[i]][['I_01_f']],
                                        I_001_m = 0,
                                        I_001_l = initial_states[[i]][['I_01_l']],
                                        
                                        C_110_u = 0,
                                        C_110_f = 0,
                                        
                                        I_110_f = 0,
                                        I_110_m = 0,
                                        I_110_l = 0,
                                        
                                        C_101_u = 0,
                                        C_101_f = 0,
                                        
                                        I_101_f = 0,
                                        I_101_m = 0,
                                        I_101_l = 0,
                                        
                                        C_011_u = initial_states[[i]][['C_11_u']],
                                        C_011_f = initial_states[[i]][['C_11_f']],
                                        
                                        I_011_f = initial_states[[i]][['I_11_f']],
                                        I_011_m = 0,
                                        I_011_l = initial_states[[i]][['I_11_l']],
                                        
                                        C_111_u = 0,
                                        C_111_f = 0,
                                        
                                        I_111_f = 0,
                                        I_111_m = 0,
                                        I_111_l = 0,
                                        
                                        D = initial_states[[i]][['D']],
                                        dTrans = initial_states[[i]][['dTrans']])
  
  
  
  
  # When intriducing as a middle drug, resistance to the current front line and last line drugs are conserved, and front line drug remains front line, last line drug remains last line
  # We basically add a zero in the middle
  # so 00 -> 000 , 10 -> 100, 01 -> 001, 11 -> 101, and all 010, 110, 011, 111 = 0
  
  initial_states_middle_new[[i]]<-    c(C_000_u = initial_states[[i]][['C_00_u']],
                                        C_000_f = initial_states[[i]][['C_00_f']],
                                        
                                        I_000_f = initial_states[[i]][['I_00_f']],
                                        I_000_m = 0,
                                        I_000_l = initial_states[[i]][['I_00_l']],
                                        
                                        C_100_u = initial_states[[i]][['C_10_u']],
                                        C_100_f = initial_states[[i]][['C_10_f']],
                                        
                                        I_100_f = initial_states[[i]][['I_10_f']],
                                        I_100_m = 0,
                                        I_100_l = initial_states[[i]][['I_10_l']],
                                        
                                        
                                        
                                        C_010_u = 0,
                                        C_010_f = 0,
                                        
                                        I_010_f = 0,
                                        I_010_m = 0,
                                        I_010_l = 0,
                                        
                                        C_001_u = initial_states[[i]][['C_01_u']],
                                        C_001_f = initial_states[[i]][['C_01_f']],
                                        
                                        I_001_f = initial_states[[i]][['I_01_f']],
                                        I_001_m = 0,
                                        I_001_l = initial_states[[i]][['I_01_l']],
                                        
                                        C_110_u = 0,
                                        C_110_f = 0,
                                        
                                        I_110_f = 0,
                                        I_110_m = 0,
                                        I_110_l = 0,
                                        
                                        C_101_u = initial_states[[i]][['C_11_u']],
                                        C_101_f = initial_states[[i]][['C_11_f']],
                                        
                                        I_101_f = initial_states[[i]][['I_11_f']],
                                        I_101_m = 0,
                                        I_101_l = initial_states[[i]][['I_11_l']],
                                        
                                        C_011_u = 0,
                                        C_011_f = 0,
                                        
                                        I_011_f = 0,
                                        I_011_m = 0,
                                        I_011_l = 0,
                                        
                                        C_111_u = 0,
                                        C_111_f = 0,
                                        
                                        I_111_f = 0,
                                        I_111_m = 0,
                                        I_111_l = 0,
                                        
                                        D = initial_states[[i]][['D']],
                                        dTrans = initial_states[[i]][['dTrans']])
  
  
  
  # When introducing drug as last line
  # All treated with former last line swithc to the new last line
  # so no one being treated with middle drug (former last line one)
  # But keep this former last line as an intermediate before escalating to the new last line one
  # 00 = 000, 10 = 100, 01 = 010, 11 = 110 (add a zero to the right of all resistance profiles)
  # 001, 011, 101 and 111 = zero because no one is resistant to this new last line drug
  
  initial_states_last_new[[i]]<-   c(C_000_u = initial_states[[i]][['C_00_u']],
                                     C_000_f = initial_states[[i]][['C_00_f']],
                                     
                                     I_000_f = initial_states[[i]][['I_00_f']],
                                     I_000_m = 0,
                                     I_000_l = initial_states[[i]][['I_00_l']],
                                     
                                     C_100_u = initial_states[[i]][['C_10_u']],
                                     C_100_f = initial_states[[i]][['C_10_f']],
                                     
                                     I_100_f = initial_states[[i]][['I_10_f']],
                                     I_100_m = 0,
                                     I_100_l = initial_states[[i]][['I_10_l']],
                                     
                                     
                                     
                                     C_010_u = initial_states[[i]][['C_01_u']],
                                     C_010_f = initial_states[[i]][['C_01_f']],
                                     
                                     I_010_f = initial_states[[i]][['I_01_f']],
                                     I_010_m = 0,
                                     I_010_l = initial_states[[i]][['I_01_l']],
                                     
                                     C_001_u = 0,
                                     C_001_f = 0,
                                     
                                     I_001_f = 0,
                                     I_001_m = 0,
                                     I_001_l = 0,
                                     
                                     C_110_u = initial_states[[i]][['C_11_u']],
                                     C_110_f = initial_states[[i]][['C_11_f']],
                                     
                                     I_110_f = initial_states[[i]][['I_11_f']],
                                     I_110_m = 0,
                                     I_110_l = initial_states[[i]][['I_11_l']],
                                     
                                     C_101_u = 0,
                                     C_101_f = 0,
                                     
                                     I_101_f = 0,
                                     I_101_m = 0,
                                     I_101_l = 0,
                                     
                                     C_011_u = 0,
                                     C_011_f = 0,
                                     
                                     I_011_f = 0,
                                     I_011_m = 0,
                                     I_011_l = 0,
                                     
                                     C_111_u = 0,
                                     C_111_f = 0,
                                     
                                     I_111_f = 0,
                                     I_111_m = 0,
                                     I_111_l = 0,
                                     
                                     D = initial_states[[i]][['D']],
                                     dTrans = initial_states[[i]][['dTrans']])
  
  


# RUN MODEL 

# results are stored in a list. Each element of the list correspond to a simulation (i simulations)
# Within each element of the list, there is one entry per model (no_int, adj_f, adj_l, diag_f, diag_l, new_f, new_l)

  
  # NO INTERVENTION: mod1, ini = initial_states
  runs[[i]]$no_int<- ode(y = initial_states[[i]], times = times, func = mod1, parms = pars[i,])
  runs[[i]]$no_int<- runs[[i]]$no_int[runs[[i]]$no_int[,'time'] %in% records,]
  
  
  
  # ADJUVANT: mod2 and mod3, ini = normal initial_states
  runs[[i]]$adj_f<-  ode(y = initial_states[[i]], times = times, func = mod2, parms = pars[i,])
  runs[[i]]$adj_f<- runs[[i]]$adj_f[runs[[i]]$adj_f[,'time'] %in% records,]
  
  
  runs[[i]]$adj_l<-  ode(y = initial_states[[i]], times = times, func = mod3, parms = pars[i,])
  runs[[i]]$adj_l<- runs[[i]]$adj_l[runs[[i]]$adj_l[,'time'] %in% records,]
  
  
  # DIAGNOSTIC: mod4 and mod5, ini = initial_states_front_diag and initial_states_last_diag
  runs[[i]]$diag_f<- ode(y = initial_states_front_diag[[i]], times = times, func = mod4, parms = pars[i,])
  runs[[i]]$diag_f<- runs[[i]]$diag_f[runs[[i]]$diag_f[,'time'] %in% records,]
  
  
  runs[[i]]$diag_l<- ode(y = initial_states_last_diag[[i]], times = times, func = mod5, parms = pars[i,])
  runs[[i]]$diag_l<- runs[[i]]$diag_l[runs[[i]]$diag_l[,'time'] %in% records,]
  
  
  # NEW DRUG: mod6, initial_states_front_new, initial_states_middle_new, initial_states_last_new
  runs[[i]]$new_f<-  ode(y = initial_states_front_new[[i]], times = times, func = mod6, parms = pars[i,])
  runs[[i]]$new_f<- runs[[i]]$new_f[runs[[i]]$new_f[,'time'] %in% records,]
  
  
  runs[[i]]$new_l<-  ode(y = initial_states_last_new[[i]], times = times, func = mod6, parms = pars[i,])
  runs[[i]]$new_l<- runs[[i]]$new_l[runs[[i]]$new_l[,'time'] %in% records,]
  
  runs[[i]]$new_m<-  ode(y = initial_states_middle_new[[i]], times = times, func = mod6, parms = pars[i,])
  runs[[i]]$new_m<- runs[[i]]$new_m[runs[[i]]$new_m[,'time'] %in% records,]
  
  
  
  # PART2 2: RECORD PUBLIC HEALTH MEASURES, PARAMETERS AND INITIAL CONDITIONS IN ONE BIG OUTPUT TABLE
  
  
  # Total PREVALENCE: integral of the total number of infected people (whatever resistance profile) in the SYMPTOMATIC (treated) class
  treated<- c('I_00_l', 'I_10_l', 'I_01_l' ,'I_11_l', 'I_00_f', 'I_10_f' ,'I_01_f', 'I_11_f')
  treated_new<- c('I_000_f', 'I_100_f', 'I_010_f', 'I_001_f', 'I_110_f', 'I_101_f', 'I_011_f', 'I_111_f',
                  'I_000_m', 'I_100_m', 'I_010_m', 'I_001_m', 'I_110_m', 'I_101_m', 'I_011_m', 'I_111_m',
                  'I_000_l', 'I_100_l', 'I_010_l', 'I_001_l', 'I_110_l', 'I_101_l', 'I_011_l', 'I_111_l')
  
  
  # Total number of DEATH: Integral of the number of deads of the course of the simulation
  deads<- 'D'  
  
  
  # No. TRANSMISSONS
  trans<- 'dTrans'   
  
  
  # PREVALENCE OF MULTIDRUG RESISTANCE
  multis<- c('C_11_u', 'C_11_f', 'I_11_f', 'I_11_l')
  multis_new<- c('C_111_u','C_111_f', 'I_111_f', 'I_111_m', 'I_111_l')
  
  
  c(a = 2,
    b = 3,
    c=4)
  
  
  out_data[i,]<- c(
    
    
    tot_prevalence_no_int = integral(runs[[i]]$no_int, treated),
    tot_prevalence_adj_f = integral(runs[[i]]$adj_f, treated),
    tot_prevalence_adj_l = integral(runs[[i]]$adj_l, treated),
    tot_prevalence_diag_f = integral(runs[[i]]$diag_f, treated),
    tot_prevalence_diag_l = integral(runs[[i]]$diag_l, treated),
    tot_prevalence_new_f = integral(runs[[i]]$new_f, treated_new),
    tot_prevalence_new_l = integral(runs[[i]]$new_l, treated_new),
    tot_prevalence_new_m = integral(runs[[i]]$new_m, treated_new),
  
  
    nb_death_no_int = integral(runs[[i]]$no_int, deads),
    nb_death_adj_f = integral(runs[[i]]$adj_f, deads),
    nb_death_adj_l = integral(runs[[i]]$adj_l, deads),
    nb_death_diag_f = integral(runs[[i]]$diag_f, deads),
    nb_death_diag_l = integral(runs[[i]]$diag_l, deads),
    nb_death_new_f = integral(runs[[i]]$new_f, deads),
    nb_death_new_l = integral(runs[[i]]$new_l, deads),
    nb_death_new_m = integral(runs[[i]]$new_m, deads),
    
    total_trans_no_int = integral(runs[[i]]$no_int, trans),
    total_trans_adj_f = integral(runs[[i]]$adj_f, trans),
    total_trans_adj_l = integral(runs[[i]]$adj_l , trans),
    total_trans_diag_f = integral(runs[[i]]$diag_f, trans),
    total_trans_diag_l = integral(runs[[i]]$diag_l, trans),
    total_trans_new_f = integral(runs[[i]]$new_f, trans),
    total_trans_new_l = integral(runs[[i]]$new_l, trans),
    total_trans_new_m = integral(runs[[i]]$new_m, trans),
    
    
    # INFECTION LENGHT: total_prevalence[i] / total_trans[i]
    
    infection_lenght_no_int = integral(runs[[i]]$no_int, treated)/integral(runs[[i]]$no_int, trans),
    infection_lenght_adj_f = integral(runs[[i]]$adj_f, treated)/integral(runs[[i]]$adj_f, trans),
    infection_lenght_adj_l = integral(runs[[i]]$adj_l, treated)/integral(runs[[i]]$adj_l, trans),
    infection_lenght_diag_f = integral(runs[[i]]$diag_f, treated)/integral(runs[[i]]$diag_f, trans),
    infection_lenght_diag_l = integral(runs[[i]]$diag_l, treated)/integral(runs[[i]]$diag_l, trans),
    infection_lenght_new_f = integral(runs[[i]]$new_f, treated_new)/integral(runs[[i]]$new_f, trans),
    infection_lenght_new_l = integral(runs[[i]]$new_l, treated_new)/integral(runs[[i]]$new_l, trans),
    infection_lenght_new_m = integral(runs[[i]]$new_m, treated_new)/integral(runs[[i]]$new_m, trans),
    
    
    # PROBA DEATH: total_death[i] / total_trans[i]
    
    proba_death_no_int = integral(runs[[i]]$no_int, deads)/integral(runs[[i]]$no_int, trans),
    proba_death_adj_f = integral(runs[[i]]$adj_f, deads)/integral(runs[[i]]$adj_f, trans),
    proba_death_adj_l = integral(runs[[i]]$adj_l, deads)/integral(runs[[i]]$adj_l, trans),
    proba_death_diag_f = integral(runs[[i]]$diag_f, deads)/integral(runs[[i]]$diag_f, trans),
    proba_death_diag_l = integral(runs[[i]]$diag_l, deads)/integral(runs[[i]]$diag_l, trans),
    proba_death_new_f = integral(runs[[i]]$new_f, deads)/integral(runs[[i]]$new_f, trans),
    proba_death_new_l = integral(runs[[i]]$new_l, deads)/integral(runs[[i]]$new_l, trans),
    proba_death_new_m = integral(runs[[i]]$new_m, deads)/integral(runs[[i]]$new_m, trans),
    
    
    multi_no_int = integral(runs[[i]]$no_int, multis),
    multi_adj_f  = integral(runs[[i]]$adj_f, multis),
    multi_adj_l = integral(runs[[i]]$adj_l, multis),
    multi_diag_f = integral(runs[[i]]$diag_f, multis), 
    multi_diag_l = integral(runs[[i]]$diag_l, multis) ,
    multi_new_f = integral(runs[[i]]$new_f, multis_new),
    multi_new_m = integral(runs[[i]]$new_m, multis_new),
    multi_new_l = integral(runs[[i]]$new_l, multis_new),
    
    
    prop_f_usage = sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f','I_00_f', 'I_10_f', 'I_01_f', 'I_11_f')])/sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f','I_00_f', 'I_10_f', 'I_01_f', 'I_11_f','I_00_l', 'I_10_l', 'I_01_l', 'I_11_l')]),
    
    prop_inapropriate_usage = sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f')])/sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f','I_00_f', 'I_10_f', 'I_01_f', 'I_11_f')]),
    
    
    
    # PARAMETERS
    
    bw = pars[i,'bw'],
    bb = pars[i,'bb'],
    c = pars[i,'c'],
    d = pars[i,'d'],
    g = pars[i,'g'],
    z = pars[i,'z'],
    a = pars[i,'a'],
    e = pars[i,'e'],
    p = pars[i,'p'],
    mu = pars[i,'mu'],
    f = pars[i,'f'],
    seed = pars[i,'my_seed'],
    draw = pars[i,'draw'],
    
    
    
    # INITIAL CONDITIONS = END OF BURN IN PHASE
    
    ini_C_00_u = initial_states[[i]][[1]],
    ini_C_00_f = initial_states[[i]][[2]],
    ini_I_00_f = initial_states[[i]][[3]],
    ini_I_00_l = initial_states[[i]][[4]],
    ini_C_10_u = initial_states[[i]][[5]],
    ini_C_10_f = initial_states[[i]][[6]],
    ini_I_10_f = initial_states[[i]][[7]],
    ini_I_10_l = initial_states[[i]][[8]],
    ini_C_01_u = initial_states[[i]][[9]],
    ini_C_01_f = initial_states[[i]][[10]],
    ini_I_01_f = initial_states[[i]][[11]],
    ini_I_01_l = initial_states[[i]][[12]],
    ini_C_11_u = initial_states[[i]][[13]],
    ini_C_11_f = initial_states[[i]][[14]],
    ini_I_11_f = initial_states[[i]][[15]],
    ini_I_11_l = initial_states[[i]][[16]],
    ini_D = initial_states[[i]][[17]],
    ini_dTrans = initial_states[[i]][[18]],
    
    
    neg_values_burn_in = sum(tail(outs_burn_in[[i]],1) < 0) > 0 ,
    neg_values_run = sum(c(tail(runs[[i]]$no_int,1) < 0, tail(runs[[i]]$adj_f,1) < 0, tail(runs[[i]]$adj_l,1) < 0,
                           tail(runs[[i]]$diag_f,1) < 0, tail(runs[[i]]$diag_l,1) < 0, tail(runs[[i]]$new_f,1) < 0,
                           tail(runs[[i]]$new_m,1) < 0, tail(runs[[i]]$new_l,1) < 0)) > 0
    
    
  )
  

  
  save.image(paste('sims_nb_', nb_sim, '_seed_', my_seed, '.RData', sep=''))
  
  write.table(out_data, 
              paste('out_data_', nb_sim, '_seed_', my_seed, '.txt',sep=''),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  print(paste('run ',i, ' done!', sep=''))
  
  i = i+1
  
}