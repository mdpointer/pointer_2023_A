# pointer_2023_A

# Data files and scripts relating to the investigation of the genetic architecture of dispersal behaviour in T.castaneum

## Contents
### Four .csv data files:
dispersal_selection_data.csv

dispersal_selection_Sexing_data.csv

F1_dispersal_assay.csv

inbred_line_dispersal_assay.csv
  
  
### Ten .R scripts
dispersal_selection_analysis.R

sex_of_dispersers.R

inbred_line_disp.R

test_dispersal_prob_equations.R

dispersal_sim_functions_nlocus.R

dispersal_sim_run.R

comparing_expt_to_sim_data.R

dispersal_sim_crosses.R

inbreeding_sim_nloci_extra_gen.R

inbreeding_sim_analysis.R
  
  
## Description of .R scripts
### dispersal_selection_analysis.R
Uses data in "dispersal_selection_data.csv" to analyse a 5 generation artificial selection experiment on dispersal propensity. Includes data wrangling, modelling, plotting and comparison of 1- and 3-day dispersal assays


### sex_of_dispersers.R
Uses data in "dispersal_selection_Sexing_data.csv" to analyse the sex ratios among dispersal phenotypes


### test_dispersal_prob_equations.R
Tests the dispersal probability equations used in the dispersal simulation model


### dispersal_sim_functions_nlocus.R
Defines an agent based simulation that models the T.castaneum artificial selection experiment


### dispersal_sim_run.R
Defines a region of parameter space, then runs the simulation using "dispersal_sim_functions_nlocus.R"


### comparing_expt_to_sim_data.R
Brings together data from "dispersal_selection_data.csv" and the .rds objects created as the outputs from "dispersal_sim_run.R" to compare experimental data to simulated scenarios


### dispersal_sim_crosses.R
Simulates crosses between simulated high and low dispersal lines


### inbreeding_sim_nloci_extra_gen.R
Simulates the inbreeding design used to generate experimentally inbred T.castaneum lines


### inbreeding_sim_analysis.R
Uses the output from "inbreeding_sim_nloci_extra_gen.R" to compare the dispersal of lines inbred under different assumes genetic architectures


### inbred_line_disp.R
Uses functions from "dispersal_sim_functions_nlocus.R" to simulate the behaviour of cross offspring of simulated dispersal lines, then compares these to the dispersal behaviour of crosses between experimental dispersal lines, using data from "inbred_line_dispersal_assay.csv"

