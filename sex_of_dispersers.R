# Looking at sex of dispersers across selection

dd <- read.csv("dispersal_selection_Sexing_data.csv")
dd <- dd %>%
  filter( disp== 0 & level=="low" |
          disp==3 & level=="high") %>% 
  rowwise %>% 
  mutate(
    measurement = ifelse( measurement == "na", NA, measurement),
    level = as.factor( level),
    sex = as.factor(sex),
    gen = as.factor(gen),
    measurement = as.integer(measurement),
    line_id = paste(line, level, sep="_"),
    level = ifelse( level=="high", "High", "Low")
  ) %>%
  drop_na(measurement) %>%
  spread(sex, measurement) %>%
  mutate( prop_m = m / (m+f),
          sample_size = m+f
          )

 # chi squared on gen 1
high_g0 <- dd %>%
  filter(gen==0,
         m!=2) %>% 
  group_by(level) %>% 
  summarise(total_f = sum(f),
            total_m = sum(m)) %>% 
  filter( level=="High") %>% 
  select(-level)

chisq.test( high_g0, p= c(0.5, 0.5))

low_g0 <- dd %>%
  filter(gen==0) %>% 
  group_by(level) %>% 
  summarise(total_f = sum(f),
            total_m = sum(m)) %>% 
  filter( level=="Low") %>% 
  select(-level)

chisq.test( low_g0, p= c(0.5, 0.5))


# binomial GLM on sex ratios in high and low regime individually

high_gall <- dd %>% 
  filter(level=="High")

glm_intercept <- glm( prop_m ~ as.numeric(gen),
                      data=high_gall, weights=sample_size, family=binomial)
summary(glm_intercept)


low_gall <- dd %>% 
  filter(level=="High")

glm_intercept <- glm( prop_m ~ as.numeric(gen),
                      data=low_gall, weights=sample_size, family=binomial)
summary(glm_intercept)

# Generate figure S2
sex_ratio_plot <- dd %>% 
  mutate(level= case_when( level=="High" ~ "High dispersal lines",
                           level=="Low" ~ "Low dispersal lines")) %>% 
  ggplot( aes( x=as.numeric(gen), y=prop_m, fill=gen)) +
  geom_hline(yintercept=0.5, alpha=0.9, linetype="dashed") +
  geom_jitter(width=0.25, alpha=0.4) +
  geom_boxplot(outlier.shape=NA, alpha=0.5) +
  scale_fill_manual(values=c("red3", "orange1","chartreuse4", "cyan4","deeppink4"))+
  theme_bw() +
  facet_wrap(~level) +
  labs(x="Generation",
       y="Proportion of males") +
  theme(legend.position = "none",
        text = element_text(size = 20))
ggsave("sex_ratio_plot.png", dpi=300)
