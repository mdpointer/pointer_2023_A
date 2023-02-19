library(tidyverse)
library(ggpubr)
library(lme4)
library(RColorBrewer)

## Experimental data
dd <- read.csv("dispersal_selection_data.csv")

dd_disp_indiv <- dd %>%
  filter(value!="na", day==4) %>% 
  group_by(gen, level, line) %>% 
  mutate(total=sum(as.numeric(value))) %>%
  group_by(gen, level, line, no_dispersals) %>%
  summarise(total=mean(total), measurement = sum(as.numeric(value))) %>%
  ungroup() %>% 
  mutate( prop = measurement / total) %>% 
  mutate( disp_x_meas = no_dispersals * measurement) %>%
  group_by(gen, level, line) %>% 
  summarise(total_disp=sum(disp_x_meas),
            total_indiv=mean(total)) %>% 
  mutate(disp_indiv = total_disp / total_indiv,
         block= ifelse(line %in% 1:8, 1, 2),
         line_id= paste(line, level, sep="_"),
         level=as.factor(level),
         line=as.factor(line),
         block=as.factor(block),
         gen = as.numeric(gen)   )


runs <- c("sum_mean_parascan_output_sexlinked_1locus.rds",
          "sum_mean_parascan_output_unsexlinked_1locus.rds",
          "sum_mean_parascan_output_unsexlinked_3locus.rds",
          "sum_mean_parascan_output_unsexlinked_5locus.rds",
          "sum_mean_parascan_output_unsexlinked_10locus.rds")


## simulated data
dd_new <- tibble(
  runs = runs) %>% 
  mutate(
  data = map(runs, ~read_rds(.x) %>% 
  filter(A!=1) %>% 
  unnest(cols=c(sum_mean_out)) %>%
  rename(line = sim,
         disp_indiv = mean_disp) %>% 
  mutate(line_id = paste(line, level, sep="_"),
         scenario = paste(h, d, A)
  )
  ) )


## COMPARATIVE MODELLING
# Fit model, gen as X^2
lmx2 <- lmer( disp_indiv ~ poly(gen,2)*level + (1 | line_id),
              data = dd_disp_indiv)


r_squared <- function(vals, preds) {
  1 - (sum((vals - preds)^2) / sum((vals - mean(preds))^2))
}

assess_a_scenario <- function(model, scenario_data) {
  pred_vals <- predict(model, newdata = scenario_data, re.form=~0)
  r2 <- r_squared(scenario_data$disp_indiv, pred_vals)
  scenario <- paste(scenario_data$h[1], scenario_data$d[1], scenario_data$A[1], sep="_")
  out <- paste(r2, scenario) #list(r2 = r2, scenario = scenario)
  return(out)
}

assess_a_dataset <- function(dataset, model) {
  dataset %>% 
    group_by(scenario) %>% 
    group_map(., ~assess_a_scenario(model, .x))
}
    
dd_new <- dd_new %>% 
  mutate(
assess_dataset_output = .$data %>% 
  map(~assess_a_dataset(.x, lmx2) )
)


r2_comp <- tibble(r2 = dd_new$assess_dataset_output %>% unlist() # changed to "dataset" here from "scenario"
) %>% 
  separate(r2, into = c("r2", "scenario"), sep = " ") %>%
  separate(scenario, c("h", "d", "A"), sep="_") %>% 
  mutate(adj_r2 = ifelse( r2<0, 0, r2),
         r2 = as.numeric(r2),
         h = as.numeric(h),
         d = as.numeric(d),
         A = as.numeric(A),
         adj_r2 = as.numeric(adj_r2)
  )

last_function <- function( dataset) {
  tibble(r2 = dataset %>% unlist()
  ) %>% 
    separate(r2, into = c("r2", "scenario"), sep = " ") %>%
    separate(scenario, c("h", "d", "A"), sep="_") %>% 
    mutate(adj_r2 = ifelse( r2<0, 0, r2),
           r2 = as.numeric(r2),
           h = as.numeric(h),
           d = as.numeric(d),
           A = as.numeric(A),
           adj_r2 = as.numeric(adj_r2)
    )
}


dd_new <- dd_new %>% 
  mutate(final_data = map(.$assess_dataset_output, ~last_function(.x))
  )

final_output <- dd_new %>% 
  select(runs, final_data) %>% 
  unnest(final_data)

final_output %>% 
  filter(runs =="sum_mean_parascan_output_unsexlinked_1locus.rds") %>% 
  filter(adj_r2 > 0.8)
#highest for 1 loc sexlinked = 0.6 0.5 0.8
#highest for 1 loc unsexlinked = 0.6 0.5 0.8
#highest for 3 loc unsexlinked = 0.6 1 0.4
#highest for 5 loc unsexlinked = 0.8 1 0.4
#highest for 10 loc unsexlinked = 1 1 0.4


### Generating parameter scan plot
final_output <- final_output %>%
  mutate( Starting_A_freq = A,
          runs = substr(runs, 26, 50))
final_output$runs <- factor(final_output$runs, levels= c("sexlinked_1locus.rds",
                                         "unsexlinked_1locus.rds",
                                         "unsexlinked_3locus.rds",
                                         "unsexlinked_5locus.rds",
                                         "unsexlinked_10locus.rds"))
final_output %>% 
  ggplot(aes(x=h, y=d, fill=adj_r2)) +
  geom_tile(colour="white") +
  scale_x_continuous(breaks=seq(0.2, 1, 0.2)) +
  scale_y_continuous(breaks=seq(0, 1, 0.25)) +
  facet_wrap(runs~Starting_A_freq, nrow=5) +
  labs( x= "Heritability",
        y = "Dominance",
        fill= (bquote(~~ R^'2'))
        ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_viridis_c(limits = c(0,1), 
                           breaks = c(0, 0.25, 0.50, 0.75, 1.00))

ggsave("parameter_scan.png", dpi=300, width=170, height=100, units="mm")




### compare sex-linked and unsexlinked 1-locus models
summarise_mean_output_by_sex <- function(input_table)
{
  input_table %>% 
    group_by(sim, gen, level, sex) %>% 
    summarise( mean_disp = mean(disp_prop) ) %>%
    ungroup() %>% 
    return(.)
}

summarise_mean_output_by_sex_for_whole_run <- function(dataframe){
out <- dataframe %>% 
  mutate(
summary = map(.$final_output, ~summarise_mean_output_by_sex(.x))
) %>% 
  select(-final_output)
return(out)
}

datasets <- c("new_parameter_scan_output_sexlinked_1locus.rds",
              "new_parameter_scan_output_unsexlinked_1locus.rds")
data_full <- tibble(
  runs = datasets,
  data = map(runs, ~read_rds(.x))
) %>% 
  mutate(
    summary_data = map(.$data, ~summarise_mean_output_by_sex_for_whole_run(.x))
  ) %>% 
  select(-data) %>% 
  mutate( unnested_summary_data = map(.$summary_data, ~unnest(.x, cols= c(summary)))) %>% 
  select(-summary_data) %>% 
  unnest(cols = c(unnested_summary_data)) %>% 
  pivot_wider(names_from = sex, values_from = mean_disp) %>% 
  filter( A != 1,
          gen ==4) %>% 
  group_by( runs,
            h, d, A,
            level) %>% 
  summarise(mean_F = mean(`F`),
            mean_M = mean(M)) %>% 
  mutate( sex_diff = mean_M - mean_F )

data_full %>%
  filter(level=="low") %>% 
  pull(sex_diff) %>% 
  max()





