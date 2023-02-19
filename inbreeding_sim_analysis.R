library(tidyverse)
library(ggpubr)
library(diptest)
rm(list = ls())

# Run the sim
#source("inbreeding_sim_nloci.R")
source("inbreeding_sim_nloci_extra_gen.R")
# 
# loci <- c(1, 3, 5, 10)
# best_h <- c(0.6, 0.6, 0.8, 1)
# best_d <- c(0.5, 1, 1, 1)
# best_A <- c(0.8, 0.4, 0.4, 0.4)
# parameter_table <- tibble(loci, best_h, best_d, best_A)
# 
# output <- parameter_table %>%
#   mutate(
#     inbreeding_sim_output = pmap(list(loci,best_h, best_d, best_A),
#                                  ~inbreeding_sim(n_loci=..1,
#                                                  h=..2,
#                                                  d=..3,
#                                                  A=..4,
#                                                  n_sim = 250,
#                                                  n_generations = 10)
#     )
#   )
# 
# 
# # manipulate data for plotting & add to output dataframe
# extra_gen <- TRUE
# n_generations <- 10
# if(extra_gen == TRUE){n_generations <- n_generations +1}
# bottleneck_gens <- c(0:9)
# n_sim <- 250

# output <- output %>%
#   mutate(
#     bar    = pmap( list(.$inbreeding_sim_output, .$loci),
#                    ~data_for_bar(..1, n_generations, ..2)), # changed this and not sure it will work
#     hist   = map( .$inbreeding_sim_output,
#                   ~data_for_hist(.x, n_generations)),
#     spider = map( .$inbreeding_sim_output,
#                   ~data_for_spider(.x))
#   )

# output %>%
#   write_rds("inbreeding_sim_output_processed")
output <- read_rds("inbreeding_sim_output_processed")

# Plots for 1-locus
plots_1locus <- ggarrange(
  plot_histogram_disp(output, 1),
  plot_barplot(output, 1),
  #plot_spiderplot(output, 1, n_generations, bottleneck_gens),
  plot_histogram_allele(output, 1),
  ncol=3, nrow=1, widths=c(3, 3, 2.5)
)

# Plots for 3-locus
plots_3locus <- ggarrange(
  plot_histogram_disp(output, 2),
  plot_barplot(output, 2),
  #plot_spiderplot(output, 2, n_generations, bottleneck_gens),
  plot_histogram_allele(output, 2),
  ncol=3, nrow=1, widths=c(3, 3, 2.5)
)

# Plots for 5-locus
plots_5locus <- ggarrange(
  plot_histogram_disp(output, 3),
  plot_barplot(output, 3),
  #plot_spiderplot(output, 3, n_generations, bottleneck_gens),
  plot_histogram_allele(output, 3),
  ncol=3, nrow=1, widths=c(3, 3, 2.5)
)

# Plots for 10-locus
plots_10locus <- ggarrange(
  plot_histogram_disp(output, 4),
  plot_barplot(output, 4),
  #plot_spiderplot(output, 4, n_generations, bottleneck_gens),
  plot_histogram_allele(output, 4),
  ncol=3, nrow=1, widths=c(3, 3, 2.5)
)

# Plot for experimental data
expt_dd <- read_csv("consec_disp_and_surf.csv") %>% 
  mutate( id = 1)
plots_experimental <- ggarrange(
  ggplot(expt_dd, aes(x=disp, fill=factor(id))) +
    geom_histogram(aes(y=..density..),alpha=0.4, position="identity", color="darkslategrey") +
    geom_density(alpha=0.4) +
    labs(x="Number dispersers",
         y = "Proportion of lines",
         fill="Gen") +
    scale_fill_manual(values=c("red", "blue")) +
    scale_x_continuous(limits=c(0, 200)) +
    scale_y_continuous(limits=c(0, 0.011), breaks=seq(0, 1, 0.002)) +
    theme_bw() +
    guides(shape = guide_legend(override.aes = list(size = 0.2)),
           color = guide_legend(override.aes = list(size = 0.2))) +
    theme(legend.position="none",
          axis.title.x = element_blank()),
  get_legend(
    plot_histogram_disp_legend(output, 4)
  ) %>% 
    as_ggplot(),
  get_legend(
    plot_barplot_legend(output, 4)
  ) %>% as_ggplot(),
  # ggplot(expt_dd) + theme_void() + theme(panel.background = element_rect(fill = 'white')),
  ncol=4, nrow=1, widths=c(3, 1.5, 1.5, 2.5)
)

ggarrange( plots_experimental,
           plots_1locus,
           plots_3locus,
           plots_5locus,
           plots_10locus,
           ncol=1, nrow=5,
           labels=list("A: Experimental data", "B:   1-locus\nstarting A freq=0.8", "C:   3-locus\nstarting A freq=0.4", "D:   5-locus\nstarting A freq=0.4", "E:   10-locus\nstarting A freq=0.4"),
           label.y= c(0.97, 1.05, 1.05, 1.05, 1.05),
           label.x= c(-0.02, 0.008, 0.008, 0.008, 0.008),
           font.label = list(size = 13, color="black"))

ggsave("inbreeding_sim_overall.png", dpi=300, width=297, height=210, units="mm")

# ggarrange( plot_histogram_allele(output, 1) + labs(title=""),
#            plot_histogram_allele(output, 2) + labs(title=""),
#            plot_histogram_allele(output, 3) + labs(title=""),
#            plot_histogram_allele(output, 4) + labs(title=""),
#            labels = c("1-locus", "3-locus", "5-locus", "10-locus")
# )
# ggsave("inbreeding_sim_allele_freq.png", dpi=200, width=200, height=130, units="mm")



#### dip.test of multimodality disperser number
library(diptest)
# 1_locus
output$hist[[1]] %>% 
  filter(gen== "gen 11") %>%
  .$n %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)
# 3_locus
output$hist[[2]] %>% 
  filter(gen== "gen 11") %>%
  .$n %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)
# 5_locus
output$hist[[3]] %>% 
  filter(gen== "gen 11") %>%
  .$n %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)
# 10_locus
output$hist[[4]] %>% 
  filter(gen== "gen 11") %>%
  .$n %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)

#### dip.test of multimodality alelle frequency
library(diptest)
# 1_locus
output$spider[[1]] %>% 
  filter(gen== 11) %>%
  pull(gen_A_freq) %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)
# 3_locus
output$spider[[2]] %>% 
  filter(gen== 11) %>%
  pull(gen_A_freq) %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)
# 5_locus
output$spider[[3]] %>% 
  filter(gen== 11) %>%
  pull(gen_A_freq) %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)
# 10_locus
output$spider[[4]] %>% 
  filter(gen== 11) %>%
  pull(gen_A_freq) %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)


# submsample 1-locus to experimental N
subsample_hist_disp <- plot_histogram_disp_sub(output, 1) + labs(title="A")

subsample_hist_freq <- plot_histogram_allele_sub(output, 1) + labs(title="B")

# ggarrange( subsample_hist_disp,
#            subsample_hist_freq)
# ggsave("inbreeding_sim_subsample.png", dpi=200, width=200, height=100, units="mm")

# subsampled sim disperser number diptest
output$hist[[1]] %>% 
  filter(sim %in% 1:64) %>% 
  filter(gen== "gen 11") %>%
  .$n %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)

# subsampled sim alelle freq diptest
output$spider[[1]] %>%
  filter(sim %in% 1:64) %>% 
  filter(gen==11) %>% 
  count(gen_A_freq) %>% 
  pull(gen_A_freq) %>% 
  dip.test(., simulate.p.value =TRUE, B=10000)




  