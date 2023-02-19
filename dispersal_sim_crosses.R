# rm(list = ls())
library(tidyverse)
library(lme4)
library(lmerTest)
# Need to have the raw simulation outputs ready to load in
# At the moment, this just takes 1 individual of each sex and each level and
# mates them in the 3 different ways - could be changed to use different ones
n_individuals_to_assay <- 200
n_sim <- 50


source("dispersal_sim_functions_1locus.R")
h = 0.6                              # trait heritability
d = 0.5                                # dominance
A = 0.8                              # starting A allele frequency
#sex_linkage = FALSE                   # "TRUE" / "FALSE"
n_sim = 50                            # number of simulations
n_generations = 4                    # generations simulated

n_inds = 200                         # population size (initial and each gen)
sex_ratio = 0.5                      # proportion of females in initial pop
low_line_selection_threshold = 0     # No emigrations required for low lines
high_line_selection_threshold = 3    # No emigrations required for high lines
ad_in = 30                           # max population size for Ad/in
pair_prod_mean = 100                 # mean of distribution of offspring prod
pair_prod_sd = 25                    # sd of distribution of offspring prod
max_mates = 3                        # maximum number of mates for any 1 female

cross_design <- tibble(
  cross_type = rep(c("Hm-Hf", "Lm-Hf", "Hm-Lf", "Lm-Lf"), each=n_sim),
  female_level = rep(c("high", "high", "low", "low"), each=n_sim),
  male_level= rep(c("high", "low", "high", "low"), each=n_sim),
  sim = rep(c(1:n_sim), 4)
)
# define functions to do crosses across architectures
read_and_filter <- function(filename){
  out <- read_rds(filename) %>% 
    filter(h >0.5 & h <0.7) %>% 
    filter(d==0.5, A==0.8) %>% 
    unnest(final_output)
  return(out)
}

get_adults_for_cross <- function(dataframe){
  out <- dataframe %>% 
    filter(gen== max(gen)) %>% 
    group_by(sim, level, sex) %>% 
    sample_n(5) %>% 
    ungroup()
  return(out)
}

perform_a_cross <- function(cross_type, female_level, male_level, sim, selection_output, linkage) {
  selection_output %>%
    filter(sim== {sim}, sex == "F", level == {female_level}) %>%
    select(genos, sex) %>%
    mutate(mates = 1) %>%
    mutate(m_genos = selection_output %>%
             filter(sim == {sim},
                    sex == "M",
                    level == {male_level}) %>% .$genos) %>%
    rowwise() %>%
    mutate(n_offspring = round(rnorm(1, pair_prod_mean, pair_prod_sd))) %>%
    mutate(n_offspring = ifelse(n_offspring < 0, 0, n_offspring)) %>%
    mutate(
      allele_1 = map2(genos, n_offspring,
                      ~ sample(
                        unlist(str_split(.x, pattern = "")), size = .y, replace = T
                      )),
      allele_2 = map2(m_genos, n_offspring,
                      ~ sample(
                        unlist(str_split(.x, pattern = "")), size = .y, replace = T
                      ))
    ) %>%
    unnest(cols = c(allele_1, allele_2)) %>%
    select(allele_1, allele_2) %>%
    mutate(sex = if (linkage) {
      ifelse(allele_2 == "Y", "M", "F")
    } else {
      sample(
        c("F", "M"),
        size = length(allele_2),
        prob = c(sex_ratio, 1 - sex_ratio),
        replace = TRUE
      )
    }) %>%
    sample_n(n_individuals_to_assay) %>%
    mutate(
      genos = paste(allele_1, allele_2, sep = ""),
      genos = ifelse(genos == "aA", "Aa", genos)
    ) %>%
    mutate( l_1 = unlist(map(.$genos, ~assign_dispersal_outcome(.x, d, h))),
            disp_prop = as.numeric(l_1)
    ) %>% 
    mutate( cross_type = {cross_type},
            sim = {sim})
}

do_the_cross_for_an_architecture <- function(cross, adults, linkage){
  out <- cross %>% 
    pmap(~perform_a_cross(..1, ..2, ..3, ..4, adults, linkage))
  return(out)
}

datasets <- c("new_parameter_scan_output_sexlinked_1locus.rds",
              "new_parameter_scan_output_unsexlinked_1locus.rds")
final_output <- tibble(
  runs = datasets,
  data = map(runs, ~read_and_filter(.x))
) %>% 
  mutate( x_linkage = ifelse( runs =="new_parameter_scan_output_sexlinked_1locus.rds",
                          TRUE,
                          FALSE))

final_output <- final_output %>% 
  mutate( adults =
            map(.$data, ~get_adults_for_cross(.x)))

final_output <- final_output %>% 
  mutate( cross_output = 
  pmap(list(.$adults, .$x_linkage),
       ~do_the_cross_for_an_architecture(cross_design, ..1, ..2))
  ) %>%
  select(-data, -adults) %>% 
  unnest(cross_output) %>% 
  unnest(cross_output)
#write_rds(final_output, "final_output_of_crosses.rds")


# Plot simulated data
final_output %>% 
  mutate(sex = ifelse( sex =="F", "female", "male"),
         runs= ifelse( runs =="new_parameter_scan_output_unsexlinked_1locus.rds", "unsex-linked", "sex-linked")
  ) %>% 
  group_by(runs, sim, sex, cross_type) %>%
  summarise(disp = sum(disp_prop)) %>%
  ggplot(aes(x = cross_type, y = disp, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_point(position = position_jitterdodge(), alpha = 0.4) +
  facet_wrap(~runs) +
  theme_bw() +
  scale_fill_manual(values = c("darkmagenta", "darkolivegreen")) +
  labs(title = " ",
       x = "Cross type",
       y = "Number of dispersers")


# read in full version of data
dd <- read.csv("/Users/mdp/Google Drive/WORK/tribolium_PhD/CH.2/R/dispersal_selection/dispersal_selection_F1_assay.csv") %>%
  mutate(cross = ifelse(cross=="HmHf", "Hm-Hf",
                        ifelse(cross=="HmLf", "Hm-Lf",
                               ifelse(cross=="LmLf", "Lm-Lf", "Lm-Hf")))) %>% 
  gather(sex, disp, male:female) %>%
  mutate( line= as.factor(line),
          block= as.factor(block),
          cross= as.factor(cross),
          sex= as.factor(sex),
          non_disp= m_idd - disp)

# run models on expt data
lmer1 <- lmer( disp ~ relevel(cross, ref="Hm-Lf") + sex + (1|line), # fixed effects
               data=dd)
summary(lmer1)
lmer2 <- lmer( disp ~ cross*sex + (1|line), # interaction
               data=dd)
summary(lmer2)

# make predictions
predict <- predict(lmer2, re.form=~(1|line), type="response")
pred <- data.frame(keyName=names(predict), value=predict, row.names=NULL) %>% 
  mutate(cross = dd$cross,
         sex= dd$sex) %>% 
  group_by(sex, cross) %>% 
  summarise(value=mean(value)) %>% 
  ungroup()

# Generate figure - sim violins with expt predictions
final_output %>% 
  mutate(sex = ifelse( sex =="F", "female", "male"),
         runs= ifelse( runs =="new_parameter_scan_output_unsexlinked_1locus.rds", "unsex-linked", "sex-linked")
  ) %>% 
  group_by(runs, sim, sex, cross_type) %>%
  summarise(disp = sum(disp_prop)) %>%
  ggplot(aes(x = cross_type, y = disp, fill = sex)) +
  geom_violin(alpha=0.4, position = position_dodge(width=0.68)) +
  geom_point(data=pred, aes(x=cross, y=value, colour=sex),
             position=position_dodge(width = 0.68),
             shape=1, size=3, stroke=1.8) +
  theme_bw() +
  facet_wrap(~runs) +
  scale_fill_manual(values = c("darkmagenta", "darkolivegreen")) +
  scale_colour_manual(values = c("darkmagenta", "darkolivegreen")) +
    labs(x = "Cross type",
         y = "Number of dispersers")

#ggsave("simulated_crosses.png", dpi=300)



## Quantify
lmer2 <- lmer( disp ~ cross*sex + (1|line),
               data=dd)
summary(lmer2)

# Fit observed data model to sim outputs
r_squared <- function(vals, preds) {
  1 - (sum((vals - preds)^2) / sum((vals - mean(preds))^2))
}

assess_a_scenario <- function(model, scenario_data) {
  pred_vals <- predict(model, newdata = scenario_data, re.form=~0)
  r2 <- r_squared(scenario_data$disp, pred_vals)
  out <- r2 #list(r2 = r2, scenario = scenario)
  return(out)
}



# Get sim data in correct format and pipe to functions to run the model trained on expt data
final_output %>% 
  group_by(runs, sex, sim, cross_type) %>% 
  summarise(disp = sum(disp_prop)) %>%
  ungroup() %>% 
  mutate( sex = ifelse(sex=="M", "male", "female"),
          cross = as.factor(cross_type),
          line = sim) %>% 
  group_by(runs) %>% 
  group_split() %>% 
  map( ~assess_a_scenario(lmer2, .x))


# Run both models afresh on the simulated data
final_output_m <- final_output %>% 
  group_by(runs, sex, sim, cross_type) %>% 
  summarise(disp = sum(disp_prop)) %>%
  ungroup() %>% 
  mutate( sex = ifelse(sex=="M", "male", "female"),
          cross = as.factor(cross_type),
          line = sim)
lmer1_s <- glm( disp ~ relevel(cross, ref="HmLf") + sex,
               data=final_output_m %>% filter(runs=="new_parameter_scan_output_sexlinked_1locus.rds"))
summary(lmer1_s)
lmer2_s <- glm( disp ~ relevel(cross, ref="HmLf")*sex,
               data=final_output_m %>% filter(runs=="new_parameter_scan_output_sexlinked_1locus.rds"))
summary(lmer2_s)

lmer1_u <- glm( disp ~ relevel(cross, ref="HmLf") + sex,
               data=final_output_m %>% filter(runs=="new_parameter_scan_output_unsexlinked_1locus.rds"))
summary(lmer1_u)
lmer2_u <- glm( disp ~ relevel(cross, ref="HmLf")*sex,
               data=final_output_m %>% filter(runs=="new_parameter_scan_output_unsexlinked_1locus.rds"))
summary(lmer2_u)

# Generate table of model outputs for simulated data
Parameter <- c( "(Intercept)", "Cross type [HmLf] (reference)",
                "Cross type [HmHf]", "Cross type [LmHf]",
                "Cross type [LmLf]", "sex", "Cross type [HmLf] x sex",
                "Cross type [LmHf] x sex", "Cross type [LmLf] x sex")
Estimate_sl <- c("34.46", "", "45.58", "30.26", "-14.94", "0.48", "", "29.32", "59.96", "30.08")
SE_sl <- c("1.37", "", "1.73", "1.73", "1.73", "1.26", "", "1.71", "1.71", "1.71")
p_sl <- c("<0.001", "", "<0.001", "<0.001", "<0.001", "0.695", "<0.001", "<0.001", "<0.001")
Estimate_us <- c("51.04", "", "28.33", "-0.66", "-30.19","-0.15", "", "-2.62", "1.80", "0.54")
SE_us <- c("0.69", "", "0.87", "0.87", "0.87", "0.61", "1.73", "1.73", "1.73")
p_us <- c("<0.001", "", "<0.001", "0.447", "<0.001", "0.807", "0.130", "0.298", "0.755")
model_summary_table <- tibble(Parameter, Estimate_sl, SE_sl, p_sl, Estimate_us, SE_us, p_us)
