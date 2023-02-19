# structure to shoot for:
# this script - creates the grid of parameters of interest
#             - runs the functions script with (source)
#             - maps the parameter grid over the main model function to create a df of dfs as output
#             - writes the output to an rds

# the main model function - simulates an initial population
#                         - runs the 'step 1 generation' for high/low lines n_gens times
#                         - binds the outputs after each other

# rm(list = ls())
# # # # # JUST NEED THESE HERE FOR TESTING
h = 0.6                              # trait heritability
d = 0.5                                # dominance
A = 0.8                              # starting A allele frequency
#sex_linkage = FALSE                   # "TRUE" / "FALSE"
n_sim = 3                            # number of simulations
n_generations = 3                    # generations simulated
n_inds = 200                         # population size (initial and each gen)
sex_ratio = 0.5                      # proportion of females in initial pop
low_line_selection_threshold = 0     # No emigrations required for low lines
high_line_selection_threshold = 3    # No emigrations required for high lines
ad_in = 30                           # max population size for Ad/in
pair_prod_mean = 100                 # mean of distribution of offspring prod
pair_prod_sd = 25                    # sd of distribution of offspring prod
max_mates = 3                        # maximum number of mates for any 1 female
n_loci = 1

#################################################
# PARAMETER SCAN SETUP
library(tidyverse)
rm(list = ls())
########## Un-comment-out the relevant functions script for the desired number of loci
source("dispersal_sim_functions_nlocus.R")  # use for unsexlinked, any number of loci
#source("emigration_sim_functions_1locus.R") # use for selinked single-locus

h <- seq(0.2, 1, 0.2)
d <- seq(0, 1, 0.25)
A <- seq(0.2, 1, 0.2)
parameter_grid <- expand_grid(h, d, A)

parameter_scan_output <- parameter_grid %>%
  mutate(final_output = pmap(
    list(h, d, A),
    ~ emigration_selection_sim(
      h = ..1,
      d = ..2,
      A = ..3,
      n_sim = 10,
      n_loci = 5,
      sex_linkage = FALSE
    )
  ))

# Summarise output/s
parameter_scan_output %>%
  mutate(sum_mean_out = map(final_output,  ~ summarise_mean_output(.))) %>%
  select(-final_output) %>%
  write_rds("sum_mean_parascan_output.rds")


parameter_scan_output %>%
  mutate(sum_prop_out = map(final_output,  ~ summarise_disp_prop_output(.))) %>%
  select(-final_output) %>%
  write_rds("sum_prop_parascan_output.rds")
