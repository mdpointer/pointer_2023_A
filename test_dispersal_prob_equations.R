### testing disp probability equations
library(tidyverse)
h <- c(1, 0.75, 0.5, 0.25, 0)
Dominance <- c(0, 0.25, 0.5, 0.75, 1) #c(-1, -0.5, 0, 0.5, 1)
expand_grid(h, Dominance) %>% 
  rowwise() %>% 
  mutate( Aa = ( 1-(0.5-((h*(Dominance - (1-Dominance)))/2))),
          aa = ( 0.5-(h/2) ),
          AA = ( 0.5+(h/2) )
  ) %>% 
  gather(geno, result, Aa:AA) %>% 
  ggplot( aes( x=h, y=result, group=geno, col=geno)) +
  geom_line(position = position_dodge(width=0.1)) +
  theme_bw() +
  theme(legend.position = c(0.8,0.25)) +
  facet_wrap(~Dominance,
             labeller = label_both) +
  labs( x= "Heritability",
        y= "Probability of dispersal",
        col="Genotype")

ggsave("disp_prob_plot_h_d.png", dpi=300, width=190, height=110, units="mm")
