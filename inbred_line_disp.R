library(tidyverse)

# read in inbred line dispersal data
dd <- read.csv("inbred_line_dispersal_assay.csv")

# view distribution
dd %>% 
  ggplot(aes(x=as.numeric(disp))) +
  geom_histogram(binwidth=10) +
  theme_bw()


# Plot disperser number histogram for inclusion in figure
expt_disp_data_histogram <- dd %>% 
  mutate(colcode = 1) %>% 
  ggplot(aes(x=as.numeric(disp))) +
  geom_histogram(aes(y=..density..),alpha=0.4, position="identity", color="darkslategrey", fill="red") +
  geom_density(alpha=0.4) +
  labs(title = "C",
       x="No dispersers",
       y = "Prop of lines",
       fill="Gen") +
  scale_x_continuous(limits=c(0, 200)) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = 0.2)),
         color = guide_legend(override.aes = list(size = 0.2))) +
  theme(legend.position=c(0.2, 0.63),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 11),
        legend.key.size =unit(0.8, "cm"))

# REQUIRES PLOTS TO HAVE BEEN MADE IN inbreeding_sim_analysis.R
ggarrange( subsample_hist_disp,
           subsample_hist_freq,
           expt_disp_data_histogram
)
ggsave("inbreeding_sim_subsample.png", dpi=200, width=200, height=130, units="mm")

# dip test for multimodality
diptest::dip.test(as.numeric(dd$disp))


