### Functions to run before sex-linked dispersal model


# BIG MODEL FUNCTION
emigration_selection_sim <- function(
                    h = 0.5,                              # trait heritability 0:1
                    d = 1,                                # trait dominance 0:1, 0 = a dominant
                    A = 0.5,                              # starting A allele frequency
                    sex_linkage = FALSE,                  # "TRUE" / "FALSE"
                    n_inds = 200,                         # population size (initial and each gen)
                    sex_ratio = 0.5,                      # proportion of females in initial pop
                    low_line_selection_threshold = 0,     # No emigrations required for low lines
                    high_line_selection_threshold = 3,    # No emigrations required for high lines
                    ad_in = 30,                           # max population size for Ad/in
                    pair_prod_mean = 100,                 # mean of distribution of offspring prod
                    pair_prod_sd = 25,                    # sd of distribution of offspring prod
                    max_mates = 3,                        # maximum number of mates for any 1 female
                    n_sim = 3,                            # number of simulations
                    n_generations = 4,                    # generations simulated
                    n_loci = 3                            # number of loci
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
    across( locus_1 : ncol(.), #
            ~ ifelse( .x =="aa", 0.5-(h/2) ,
                           ifelse( .x =="AA", 0.5+(h/2) ,
                                   1-(0.5-((h*(d - (1-d)))/2)) # remaining option is Aa
                           ) ) ) )
 
 output <- output %>% 
   mutate(
     prob = rowMeans(select(output, starts_with("locus"))),
     gen = current_gen,
     level = "low",
     sex = sample( c("F", "M"),
                 size= length(genos),
                 prob= c(sex_ratio, 1-sex_ratio),
                 replace=TRUE )) %>% 
   mutate(   l_1 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
             l_3 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
             l_2 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
             disp_prop = as.double( l_1 + l_2 +l_3)
             ) %>%
   select( - c(starts_with("locus")))
 
 output <- output %>% 
   mutate( level = "high") %>% 
   rbind( output, . )
    
    
    
    for (j in 1:n_generations ) {# LOW
      current_gen <- ifelse(j==1,
                            1,
                            current_gen +1
      )
      
      test_cond <- output %>%
        filter(gen== current_gen-1,
               level== "low",
               disp_prop== low_line_selection_threshold
        )
      if( nrow(test_cond)== 0) {
        break
      } else if(
        !all( c("M" %in% test_cond$sex,
                "F" %in% test_cond$sex)) ) {
        break
      }
      


      output <- step_a_generation(
        starting_individuals = output,
        high_or_low_step = "low",
        current_generation = current_gen,
        h = h,
        d = d,
        n_inds = n_inds,
        sex_linkage = sex_linkage,
        sex_ratio = sex_ratio,
        low_line_selection_threshold = low_line_selection_threshold,
        high_line_selection_threshold = high_line_selection_threshold,
        ad_in = ad_in,
        pair_prod_mean = pair_prod_mean,
        pair_prod_sd = pair_prod_sd,
        max_mates = max_mates,
        n_loci = n_loci)
    }

    for (k in 1:n_generations ) {# HIGH
      current_gen <- ifelse(k==1,
                            1,
                            current_gen +1
      )
      
      test_cond <- output %>%
        filter(gen== current_gen-1,
               level== "high",
               disp_prop== high_line_selection_threshold
        )
      if( nrow(test_cond)== 0) {
        break
      } else if(
        !all( c("M" %in% test_cond$sex,
                "F" %in% test_cond$sex)) ) {
        break
      }
      
      
      output <- step_a_generation(
                        starting_individuals = output,
                        high_or_low_step = "high",
                        current_generation = current_gen,
                        h = h,
                        d = d,
                        n_inds = n_inds,
                        sex_linkage = sex_linkage,
                        sex_ratio = sex_ratio,
                        low_line_selection_threshold = low_line_selection_threshold,
                        high_line_selection_threshold = high_line_selection_threshold,
                        ad_in = ad_in,
                        pair_prod_mean = pair_prod_mean,
                        pair_prod_sd = pair_prod_sd,
                        max_mates = max_mates,
                        n_loci = n_loci)
    }
    
    
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
    select(sex, gen, level, genos, disp_prop, sim)
  return(final_output)
}


######
#####
####
# Runs one generation of a specified 'level' - high / low
step_a_generation <- function(starting_individuals,
                              high_or_low_step,
                              current_generation,
                              h,
                              d,
                              n_inds,
                              sex_linkage,
                              sex_ratio,
                              low_line_selection_threshold,
                              high_line_selection_threshold,
                              ad_in,
                              pair_prod_mean,
                              pair_prod_sd,
                              max_mates,
                              n_loci
                              )
# starting_individuals <- output
# high_or_low_step <- "low"
# current_generation <- current_gen <- 1
{
disp_level <- ifelse (high_or_low_step == "low",
                      low_line_selection_threshold,
                      high_line_selection_threshold)
  
adults_in <-  starting_individuals %>%
  filter(gen== current_generation -1,
         level== high_or_low_step,
         disp_prop == disp_level) %>%
  sample_n( size = calc_sample_size(., disp_level, current_generation-1, ad_in))

names <- 1:n_loci %>% paste0("ofsp_", .) # these 2 lines allow passing of variable number
vars <- rlang::syms(names)               # of ofsp_ variables to function in next pipe
 
next_gen <- adults_in %>% 
    filter(sex=="F") %>% 
    select( genos, sex) %>% 
    mutate( mates = sample( c(1:max_mates),
                            size=nrow(.),
                            replace=TRUE)
            ) %>%
    rowwise() %>%
    mutate(n_offspring = round(rnorm(1, pair_prod_mean, pair_prod_sd))) %>%
    mutate(n_offspring = ifelse(n_offspring <0, 0, n_offspring)
    ) %>%
    ungroup() %>% 
    mutate( m_genos =
              map(.$mates, ~get_male_genos(.x, filter(adults_in, sex=="M") %>% .$genos))
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
        across(5:ncol(.),
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
  mutate( genome = str_replace_all(genome, "aA", "Aa")) %>% 
  mutate( locus = str_split(.$genome, pattern = " ")) %>%
  unnest_wider(locus, names_sep="_") %>% 
  mutate(
    across( locus_1 : ncol(.), #
            ~ ifelse( .x =="aa", 0.5-(h/2) ,
                           ifelse( .x =="AA", 0.5+(h/2) ,
                                   1-(0.5-((h*(d - (1-d)))/2)) # remaining option is Aa
                           ) ) ) )
 
  next_gen <- next_gen %>% 
    mutate(
  prob = rowMeans(select(., starts_with("locus_"))),
  gen = current_generation,
  level = high_or_low_step,
  sex = sample( c("F", "M"),
                size= length(prob),
                prob= c(sex_ratio, 1-sex_ratio),
                replace=TRUE )) %>% 
  sample_n(size = ifelse(nrow(.) < n_inds,
                         nrow(.),
                         n_inds )
  ) %>% 
  mutate(   l_1 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
            l_3 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
            l_2 = map_lgl(.$prob, ~assign_dispersal_outcome(.x)),
            disp_prop = as.double( l_1 + l_2 +l_3)
  ) %>%
  select( - c(starts_with("locus"))) %>% 
  rename( genos = genome)

output <- starting_individuals %>%
    rbind(next_gen)
  

  
  return(output)
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

# FUNCTION - paste together all the genotypes of an offspring into a haplotype
# create_indiv_genome <- function(a1, a2, a3){ #ofsp_1, ofsp_2, ofsp_3
#   out <- pmap(list(a1, a2, a3), #
#               ~paste(..1, ..2, ..3, #
#                      sep=" ")
#   )
#   return(out)
# }
create_indiv_genome <- function(...){ 
  out <- pmap(list(...), 
              ~paste(..., 
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
