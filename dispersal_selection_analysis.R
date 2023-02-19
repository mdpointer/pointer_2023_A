### Dispersal selection analysis
library(tidyverse)
library(ggpubr)
library(lme4)
library(modelr)
library(performance)
library(lmerTest)

# read in data
dd <- read.csv("dispersal_selection_data.csv")

## INDIVIDUAL-LEVEL DATA
# Wrangle data to get 1 row per individual, day 4 data
dd_day_4_row_per <- dd %>%
  filter(value!="na", day==4) %>% 
  group_by(gen, level, line) %>% 
  mutate(total=sum(as.numeric(value))) %>%
  group_by(gen, level, line, no_dispersals) %>%
  summarise(total=mean(total), measurement = sum(as.numeric(value))) %>%
  ungroup() %>% 
  mutate( prop = measurement / total) %>% 
  mutate( disp_x_meas = no_dispersals * measurement) %>% 
  slice(rep(1:nrow(.),.$measurement)) %>%
  select(-measurement, -total, -prop, -disp_x_meas) %>% 
  mutate( block = ifelse(line %in% 1:8, 1, 2 ) )

## MEAN DATA
# Wrangle data to get mean emigrations per individual per line
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

## MODELLING
# Fit model, gen as X^2
lmx2 <- lmer( disp_indiv ~ poly(gen,2)*level + block + (1 | line_id),
              data = dd_disp_indiv)
summary(lmx2)

#
lmx2_high <- lmer( disp_indiv ~ poly(gen, 2) + block + (1 | line_id),
              data = dd_disp_indiv %>% filter(level=="high"))
summary(lmx2_high)


# Check residuals
dd_disp_indiv %>% 
  gather_residuals(lmx2) %>% 
  ggplot(aes(gen, resid, colour=level)) +
  geom_point() +
  theme_bw() +
  facet_wrap(model~level)


## GENERATE TABLE 1
# summary table from model output (lmx2)
parameter <- c( "(Intecept)", "Generation", "Generation^2",
                "Selection regime", "block", "Generation x Selection regime",
                "Generation^2 x Selection regime")
Estimate <- c( 2.15, 2.21, -0.75, -1.10, -0.06, -5.64, 6.26)
SE <- c( 0.03, 0.30, 0.30, 0.04, 0.04, 0.42, 0.42)
p <- c( "<0.001", "<0.001", "0.013", "<0.001", "0.146", "<0.001", "<0.001")
disp_selection_table <- tibble( parameter, Estimate, SE, p)


### Analysis of 1- vs 3-opportunity data
# models on 3-opportunity data
glm_3op <- glm( disp_indiv ~ poly(gen,2) + level + block,
               data = dd_disp_indiv_4)
summary(glm_3op)
glm_3op <- glm( disp_indiv ~ poly(gen,2)*level + block,
                data = dd_disp_indiv_4)
summary(glm_3op)

# models on 1-opportunity data
glm_1op <- glm( disp_indiv ~ poly(gen,2) + level + block,
               data = dd_disp_indiv_2)
summary(glm_1op)





# 3-opportunities
dd_day_4_row_per <- dd %>%
  filter(value!="na", day==4) %>% 
  group_by(gen, level, line) %>% 
  mutate(total=sum(as.numeric(value))) %>%
  group_by(gen, level, line, no_dispersals) %>%
  summarise(total=mean(total), measurement = sum(as.numeric(value))) %>%
  ungroup() %>% 
  mutate( prop = measurement / total) %>% 
  mutate( disp_x_meas = no_dispersals * measurement) %>% 
  slice(rep(1:nrow(.),.$measurement)) %>%
  select(-measurement, -total, -prop, -disp_x_meas) %>% 
  mutate( block = ifelse(line %in% 1:8, 1, 2 ) )


dd_disp_indiv_4 <- dd %>%
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


## GENERATE FIGURE 2
# Plot means per line with model predictions
pred_plot <- dd_disp_indiv %>%
  filter(gen<6) %>% 
  ggplot( aes( x=gen+1, y=disp_indiv, col=level )) +
  geom_line( aes(group=line_id), alpha=0.4) +
  geom_point( alpha=0.4) +
  geom_smooth(method = "lm", formula=y~poly(x,2), level=0.95, linetype="dashed") +
  geom_point(data=filter(dd_disp_indiv, gen==0, level=="high"),
             aes(x=gen+1, y=disp_indiv), col="black") +
  theme_bw() +
  scale_color_manual(values=c("#D35400", "#0F4880")) +
  scale_y_continuous(limits = c(0, 3)) +
  labs(x="Generation", y= expression(Dispersals ~individual^-1),
       col = "Selection\nregime") +
  theme(plot.title = element_text(face = "bold"),
        text = element_text(size = 20))
pred_plot
#ggsave("disp_selection_pred_plot.png",dpi = 300)



# save data as .rds for use in comparison with simulated data
dd_disp_indiv %>%
  group_by(gen, level) %>% 
  summarise( mean = mean(disp_indiv),
             sd = sd(disp_indiv)) %>% 
  write_rds("experimental_values_for_comparison.rds")




# Looking at 3- vs 1-opportunity assay
dd_disp_indiv_gen4_day4 <- dd %>%
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
         gen = as.numeric(gen)   ) %>% 
  filter(gen==4)

# Fit model to 3-opportunity data
day4_model <- glm( disp_indiv ~ level + block,
              data = dd_disp_indiv_gen4_day4)
summary(day4_model)



# Looking at 1-opportunity assay
dd_disp_indiv_gen4_day2 <- dd %>%
  filter(value!="na", day==2) %>% 
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
         gen = as.numeric(gen)   ) %>% 
  filter(gen==4)

# Fit model to 1-opportunity data
day2_model <- glm( disp_indiv ~ level + block,
                   data = dd_disp_indiv_gen4_day2)
summary(day2_model)

  


  
  