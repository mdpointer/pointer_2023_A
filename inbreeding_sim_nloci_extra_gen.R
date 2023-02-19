### Functions to run before sex-linked dispersal model

# BIG MODEL FUNCTION
inbreeding_sim <- function(
  h = 0.5,                              # trait heritability 0:1
  d = 1,                                # trait dominance 0:1, 0 = a dominant
  A = 0.5,                              # starting A allele frequency
  n_inds = 200,                         # population size (initial and each gen)
  ad_in = 30,                           # max population size for Ad/in
  pair_prod_mean = 100,                 # mean of distribution of offspring prod
  pair_prod_sd = 25,                    # sd of distribution of offspring prod
  n_sim = 3,                            # number of simulations
  n_generations = 10,                   # generations simulated
  n_loci = 1,                            # number of loci
  max_mates = 3
) {
  
  A_freq <- A
  a1 <- c( "A", "a")
  a1_prob <- c(A_freq, 1 - A_freq)
  a2<- c("A", "a")
  a2_prob <- c(A_freq, 1-A_freq)
  
  for ( i in 1:n_sim) {
    #simulate initial population
    current_gen <- 0
    
    output <- tibble(locus_1 = paste(
      sample(a1,
             n_inds,
             prob = a1_prob,
             replace = T),
      sample (a2,
              n_inds,
              prob = a2_prob,
              replace = T),
      sep = ""
    ) %>% 
      str_replace("aA", "Aa" ))
    
    if( n_loci != 1) {
      for( z in 2:n_loci) { #
        new <- paste(
          sample(a1,
                 n_inds,
                 prob = a1_prob,
                 replace = T),
          sample (a2,
                  n_inds,
                  prob = a2_prob,
                  replace = T),
          sep = ""
        ) %>% str_replace("aA", "Aa" )
        output[,ncol(output) +1] <- new
        colnames(output)[ncol(output)] <- paste0("locus_", z)
      }
    }
    
    output <- output %>% unite(genos, remove=F, sep= " ") %>% 
      mutate(
        gen = current_gen) %>% 
      select(-c(starts_with("locus_")))
    
    for (j in 1:n_generations ) {
      
      current_gen <- ifelse(j==1,
                            1,
                            current_gen +1 )
      
      output <- step_a_generation(
        starting_individuals = output,
        current_generation = current_gen,
        h = h,
        d = d,
        pair_prod_mean = pair_prod_mean,
        pair_prod_sd = pair_prod_sd,
        n_loci = n_loci)
      
    }
    
    gen <- n_generations +1
    adults_in <- output %>% 
      filter(gen == n_generations)
    
    names <- 1:n_loci %>% paste0("ofsp_", .) # these 2 lines allow passing of variable number
    vars <- rlang::syms(names) 
    
    next_gen <- adults_in %>%
      sample_n(ad_in) %>% 
      select( genos ) %>% 
      mutate( mates = sample( c(1:max_mates),
                              size=nrow(.),
                              replace=TRUE)
      ) %>%
      rowwise() %>%
      mutate(n_offspring = 100) %>%
      ungroup() %>% 
      mutate( m_genos =
                map(.$mates, ~get_male_genos(.x, adults_in %>% .$genos))
      ) %>%
      unnest_wider(m_genos, names_sep="_") %>%
      rowwise() %>% 
      mutate( m_genos_3 = ifelse( !("m_genos_3" %in% colnames(.)), NA, m_genos_3),
              m_genos_2 = ifelse( !("m_genos_2" %in% colnames(.)), NA, m_genos_2)) %>% 
      replace(is.na(.),"") %>%
      mutate(
        genos = str_split(genos, " ")
      ) %>% 
      mutate(
        across(4:ncol(.),
               ~str_split(.x, " ")
        ) ) %>%
      mutate( males= list(mating_combine_males(m_genos_1, m_genos_2, m_genos_3))) %>%
      select( -m_genos_1, -m_genos_2, -m_genos_3 ) %>% 
      rowwise() %>% 
      mutate( ofsp = list(mating_combine_f_m(genos, males, n_offspring))) %>% 
      unnest_wider(ofsp, names_sep="_") %>% 
      rowwise() %>%
      mutate( genome = list(create_indiv_genome(!!!vars))) %>% 
      unnest(genome) %>% 
      select(genome) %>%
      mutate( genos = str_replace_all(genome, "aA", "Aa"),
              gen = gen) %>% 
      select(genos, gen) %>% 
      sample_n(200)
    
    output <- output %>% 
      rbind(., next_gen)
    
    output_sim <- output %>% 
      mutate(
        sim= i
      )
    
    ifelse(i==1,
           final_output <- output_sim,
           final_output <- rbind(final_output, output_sim)
    )
    
  }
  
  final_output <- final_output %>% 
    dispersal_assay_on_inbred_lines(h=h, d=d)
  
  return(final_output)
}


###### INTERIOR FUNCTIONS
######
# Runs one generation of inbreeding
step_a_generation <- function(starting_individuals,
                              current_generation,
                              h,
                              d,
                              n_inds,
                              sex_ratio,
                              ad_in,
                              pair_prod_mean,
                              pair_prod_sd,
                              n_loci
)
  # starting_individuals <- output
  # current_generation <- current_gen <- 1
{
  
  adults_in <-  starting_individuals %>%
    filter(gen== current_generation -1)
  
  names <- 1:n_loci %>% paste0("ofsp_", .) # these 2 lines allow passing of variable number
  vars <- rlang::syms(names)               # of ofsp_ variables to function in next pipe
  
  next_gen <- adults_in %>% 
    sample_n(1) %>%
    select( genos ) %>% 
    mutate( mates = 1
    ) %>%
    rowwise() %>%
    mutate(n_offspring = 200 ) %>%
    ungroup() %>% 
    mutate( m_genos =
              map(.$mates, ~get_male_genos(.x, sample_n(adults_in, 1) %>% .$genos))
    ) %>%
    mutate(
      genos = str_split(genos, " "),
      males = str_split(m_genos, " ")
    ) %>% 
    rowwise() %>% 
    mutate( ofsp = list(mating_combine_f_m(genos, males, n_offspring))) %>% 
    unnest_wider(ofsp, names_sep="_") %>% 
    rowwise() %>% 
    mutate( genome = list(create_indiv_genome(!!!vars))) %>% 
    unnest(genome) %>% 
    select(genome) %>%
    mutate( genos = str_replace_all(genome, "aA", "Aa"),
            gen = current_generation) %>% 
    select(-genome)
  
  output <- starting_individuals %>%
    rbind(next_gen)
  
  return(output)
}

##### ASSIGN DISPERSAL OUTCOMES
dispersal_assay_on_inbred_lines <- function(dataframe, h, d) {
  dispersal_outcomes <- dataframe %>% 
    mutate( locus = str_split(.$genos, pattern = " ")) %>%
    unnest_wider(locus, names_sep="_") %>%
    mutate(
      across( locus_1 : ncol(.), #
              ~ ifelse( .x =="aa", 0.5-(h/2) ,
                        ifelse( .x =="AA", 0.5+(h/2) ,
                                1-(0.5-((h*(d - (1-d)))/2)) # remaining option is Aa
                        ) ) ) ) %>%
    mutate(  prob = rowMeans(select(., starts_with("locus"))) ) %>% 
    mutate(   l_1 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
              l_3 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
              l_2 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
              disp_prop = as.double( l_1 + l_2 +l_3)
    ) %>%
    select(-c(prob, l_3, l_2, starts_with("locus_"))) %>% 
    return(dispersal_outcomes)
}



# FUNCTION - paste together male genotypes for all mates of a given female
mating_combine_males <- function(m1, m2, m3){
  out <- pmap(list(m1, m2, m3), ~paste0(..1, ..2, ..3))
  return( out)
}

# FUNCTION - sample female and male/s alleles to create offspring genotypes
mating_combine_f_m <- function(f, m, n_offsp) {
  out <- pmap(list(f, m, n_offsp), ~paste0(
    sample(
      unlist(str_split(..1, pattern="")),
      ..3,
      replace = T
    ),
    sample(
      unlist(str_split(..2, pattern="")),
      ..3,
      replace = T
    )
  )
  )
}

create_indiv_genome <- function(...){ #
  out <- pmap(list(...), #
              ~paste(..., #
                     sep=" ")
  )
  return(out)
}

#FUNCTION
# Generate dispersal outcomes based on genotype, dominance and heritability
assign_dispersal_outcome <- function(probability) {
  out <- rbernoulli(1, p = probability)
  return(out)
}

#FUNCTION
# assign sample_size for adults in, based on max possible from assay outcome
calc_sample_size <-
  function(input_table,
           filter_level,
           current_generation,
           adults_in_number) {
    if (nrow(
      filter(input_table, disp_prop ==filter_level & gen ==current_generation)
    ) >= adults_in_number) {
      return(adults_in_number)
    } else{
      return(nrow(filter(
        input_table, disp_prop == filter_level
      )))
    }
  }

#FUNCTION
# for each female, sample the male genotypes 'mates' number of times
get_male_genos <- function(times, genos) {
  sample( (genos),
          size=times,
          replace=TRUE)
}

### FUNCTIONS FOR SUMMARY AND ANALYSIS

# FUNCTION - Summarise output as means
summarise_mean_output <- function(input_table)
{
  input_table %>% 
    group_by(sim, gen, level) %>% 
    summarise( mean_disp = mean(disp_prop) ) %>%
    ungroup() %>% 
    return(.)
}

# FUNCTION - Summarise output as proportions of dispersers
summarise_disp_prop_output <- function(input_table)
{
  input_table %>% 
    group_by(sim, gen, level, disp_prop) %>% 
    summarise( value=n()) %>% 
    group_by(gen, level, disp_prop) %>% 
    summarise( value= mean(value)) %>%
    ungroup() %>%
    mutate(gen = as.numeric(gen) ) %>%
    return(.)
}

# Modified from somebody else's code for multiple histograms on same plot
plot_multi_histogram <- function(df, feature, label_column, rownumber) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(aes(y=..density..),alpha=0.4, position="identity", color="darkslategrey") +
    geom_density(alpha=0.4) +
    labs(x="Number of dispersers",
         y = "Proportion of lines",
         fill="Gen") +
    scale_fill_manual(values=c('blue', 'red')) +
    scale_x_continuous(limits=c(0, 200)) +
    theme_bw() +
    guides(shape = guide_legend(override.aes = list(size = 0.2)),
           color = guide_legend(override.aes = list(size = 0.2))) +
    if(rownumber==1){theme(legend.position="none",
                           axis.title.x = element_blank())}
  else if(rownumber==2){theme(legend.position="none",
                              axis.title.x = element_blank())}
  else if(rownumber==3){theme(legend.position="none",
                              axis.title.x = element_blank())}
  else if(rownumber==4){theme(legend.position="none")}
  plt + guides(fill=guide_legend(title=label_column))
}

### MANIPULATE DATA FOR PLOTTING
# data for barplot
data_for_bar <- function(dataframe, gens, loci) {
  nos <- c(1:loci)
  out <- dataframe %>%
    filter(gen== 11) %>%
    mutate( genos = strsplit(genos, " ")) %>%
  rowwise() %>% 
  mutate( genos = list(pmap(list(genos, nos), ~ paste(..1, ..2, sep="_")))) %>% 
  unnest_longer(genos) %>% 
  mutate( genos = str_split(genos, "_")) %>% 
  unnest_wider(genos) %>% 
    rename( genos = 1,
            locus = 2) %>% 
      group_by(sim, locus) %>% 
    distinct(genos) %>% 
    summarise(genos = paste0(genos, collapse = "")) %>%
    mutate( outcome = ifelse(str_detect(genos, "A") & str_detect(genos, "a"), "Both present",
                             ifelse(str_detect(genos, "A"), "A fixed",
                                    "a fixed"))) %>% 
    group_by(locus) %>% 
    count(outcome) %>% 
    right_join(.,
               expand.grid(locus = as.character( c(1:10) ),
                           outcome = c("a fixed", "A fixed", "Both present") ) ) %>% 
    mutate( prop= (n/n_sim)*100,
            prop = replace_na(prop, 0) )
  return(out)
}

# data for spider plot
data_for_spider <- function(dataframe){
  out <- dataframe %>%
    mutate( genos = strsplit(genos, " ")) %>% 
    unnest_longer(genos) %>% 
    group_by(sim, gen) %>%
    count(genos) %>%
    spread(key=genos, value=n) %>%
    replace_na(list(gen=0, aa=0, Aa=0, AA=0)) %>%
    mutate(gen_A_freq = (Aa + (AA*2)) / ((aa*2) + (Aa*2) + (AA*2)) )
  return(out)
}

# data for histogram
data_for_hist <- function(dataframe, gens){
  out <- dataframe %>% 
    filter(gen==gens | gen==0) %>% 
    group_by(sim, gen) %>% 
    count(l_1) %>% # change l_1 to disp_prop here if you want to look at 3-opportunity to dispersal
    filter(l_1==TRUE) %>% 
    mutate( gen = ifelse( gen==0, "gen 0", "gen 11"))
  return(out)
}

# Make plots
# Plot histogram plot
plot_histogram_disp <- function(dataframe, rownumber) {
  plot <- dataframe$hist[[rownumber]] %>% 
    plot_multi_histogram("n", "gen", rownumber)
  plot
}

# # Plot spider plot
# plot_spiderplot <- function(dataframe, rownumber, gens, bngens){
#   plot <- dataframe$spider[[rownumber]] %>% 
#     ggplot( aes( x= gen, y=gen_A_freq, group=sim)) +
#     geom_line(alpha=0.5) +
#     scale_colour_manual(values=c("firebrick", "black")) +
#     scale_x_continuous(breaks=seq(0,gens, by=2)) +
#     scale_y_continuous(limits =c(0, 1)) +
#     geom_rect(data=tibble( xmax = c(bngens +0.4),
#                            xmin= c(bngens -0.4),
#                            ymax=Inf,
#                            ymin=-Inf),
#               aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
#               fill="firebrick",
#               alpha=0.2,
#               inherit.aes = FALSE) +
#     theme_bw() +
#     labs(colour = "Generation",
#          x="Generation",
#          y='A allele frq') +
#     theme(legend.position="left") +
#     if(rownumber!=4){theme(axis.title.x = element_blank())}
#   plot
# }


# Plot barplot
plot_barplot <- function(dataframe, rownumber){
  plot <- barplot <- dataframe$bar[[rownumber]] %>% 
    ggplot( aes(x=factor(locus, levels= c(1:10)), y=prop, fill=outcome)) +
    geom_bar(alpha=0.8, stat='identity', position='dodge') +
    theme_bw() +
    scale_fill_manual(values=c("blue4", "olivedrab4", "gold3", "darkmagenta")) +
    scale_y_continuous(limits = c(0, 80)) +
    labs(x='Locus',
         y='% of runs',
         fill = "Alleles in gen 11") +
    guides(shape = guide_legend(override.aes = list(size = 0.2)),
           color = guide_legend(override.aes = list(size = 0.2))) +
    if(rownumber==1){theme(legend.position="none",
                           axis.title.x = element_blank())}
  else if(rownumber==2){theme(legend.position="none",
                              axis.title.x = element_blank())}
  else if(rownumber==3){theme(legend.position="none",
                              axis.title.x = element_blank())}
  else if(rownumber==4){theme(legend.position="none")}
  plot
}

# Plot histogram of allele F
plot_histogram_allele <- function(dataframe, rownumber) {
  plot <- dataframe$spider[[rownumber]] %>%
    filter(gen==11) %>%
    ggplot(aes(x=gen_A_freq)) +
    geom_histogram() +
    scale_y_continuous(limits=c(0,190)) +
    labs(x="A allele freq") +
    theme_bw() +
    if(rownumber!=4){theme(axis.title.x = element_blank())}
  plot
}



##Plot functions for subsampled data (N=64)
# Plot histogram of disp outcomes of subsample
plot_histogram_disp_sub <- function(dataframe, rownumber) {
  plot <- dataframe$hist[[rownumber]] %>% 
    filter(sim %in% 1:64) %>% 
    plot_multi_histogram_subsample("n", "gen", rownumber)
  plot
}

# # Plot histogram of allele F of subsample
# plot_histogram_allele_sub <- function(dataframe, rownumber) {
#   plot <- dataframe$spider[[rownumber]] %>%
#     filter(sim %in% 1:64) %>% 
#     filter(gen==11) %>%
#     ggplot(aes(x=gen_A_freq)) +
#     geom_histogram() +
#     scale_y_continuous(limits=c(0,60)) +
#     labs(x="'A' allele freq") +
#     theme_bw()
#   plot
# }  ## delete if subsample plot cut from paper

# # Modified from somebody else's code for multiple histograms on same plot
# plot_multi_histogram_subsample <- function(df, feature, label_column, rownumber) {
#   plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
#     geom_histogram(aes(y=..density..),alpha=0.4, position="identity", color="darkslategrey") +
#     geom_density(alpha=0.4) +
#     labs(x="No dispersers",
#          y = "Prop of lines",
#          fill="Gen") +
#     scale_fill_manual(values=c('blue', 'red')) +
#     scale_x_continuous(limits=c(0, 200)) +
#     theme_bw() +
#     guides(shape = guide_legend(override.aes = list(size = 0.2)),
#            color = guide_legend(override.aes = list(size = 0.2))) +
#     theme(legend.position=c(0.2, 0.63),
#                            legend.title = element_text(size = 14), 
#                            legend.text = element_text(size = 11),
#                            legend.key.size =unit(0.8, "cm"))
#   plt + guides(fill=guide_legend(title=label_column))
# }  ## delete if subsample plot cut from paper

# Extracting the legend for the histogram
plot_multi_histogram_legend <- function(df, feature, label_column, rownumber) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(aes(y=..density..),alpha=0.4, position="identity", color="darkslategrey") +
    geom_density(alpha=0.4) +
    scale_fill_manual(values=c('blue', 'red')) +
    scale_x_continuous(limits=c(0, 200)) +
    theme_bw() +
    labs(x="Number of dispersers",
         y = "Proportion of lines") +
    guides(shape = guide_legend(override.aes = list(size = 0.2)),
           color = guide_legend(override.aes = list(size = 0.2)))
  plt + guides(fill=guide_legend(title="Generation"))
}

plot_histogram_disp_legend <- function(dataframe, rownumber) {
  plot <- dataframe$hist[[rownumber]] %>% 
    plot_multi_histogram_legend("n", "gen", rownumber)
  plot
}


# Extracting the legend for the barplot
plot_barplot_legend <- function(dataframe, rownumber){
  plot <- barplot <- dataframe$bar[[rownumber]] %>% 
    ggplot( aes(x=factor(locus, levels= c(1:10)), y=prop, fill=outcome)) +
    geom_bar(alpha=0.8, stat='identity', position='dodge') +
    theme_bw() +
    scale_fill_manual(values=c("blue4", "olivedrab4", "gold3", "darkmagenta")) +
    scale_y_continuous(limits = c(0, 80)) +
    labs(x='Locus',
         y='% of runs',
         fill = "Alleles in generation 11") +
    guides(shape = guide_legend(override.aes = list(size = 0.2)),
           color = guide_legend(override.aes = list(size = 0.2)))
  plot
}




