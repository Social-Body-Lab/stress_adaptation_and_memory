"
Authors: Arran J. Davis
Emails: arran.davis@anthro.ox.ac.uk | davis.arran@gmail.com
Affiliation: Social Body Lab, Institute of Human Sciences, University of Oxford
Date: 01/08/2022
"

library(extrafont)
library(faux)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(tidyr)
library(broom.mixed)
library(purrr)
library(censReg)
library(plm)
library(reshape2)
library(data.table)

#clean environment
rm(list = ls())

#set working directory
setwd(getSrcDirectory()[1])
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load plot themes
source("plot_theme.R")

################################################################################################################################################

### CREATE DATA ###

#create the neighbourhood variable 
between = list(neighbourhood = c(Scampia = "Scampia", 
                                 Roma = "Roma camp",
                                 Pozzuoli = "Pozzuoli "))

#create the test variable
within = list(test = c("Verbal short-term memory",
                       "Visuospatial short-term memory",
                       "Verbal working memory", 
                       "Visuospatial working memory"))

within_short = list(test = c("verbstm",
                             "visstm",
                             "verbwm", 
                             "viswm"))

#create scores for each neighbourhood for verbal short-term memory
scampia_verb_stm = 4.5
roma_verb_vstm = 4.5
pozzuoli_verb_vstm = 5

#create scores for each neighbourhood for visuospatial short-term memory 
scampia_vis_stm = 4.5
roma_vis_vstm = 4.5
pozzuoli_vis_vstm = 4.8

#create scores for each neighbourhood for verbal working memory
scampia_verb_wm = 3.5
roma_verb_wm = 3.5
pozzuoli_verb_wm = 4

#create scores for each neighbourhood for visuospatial working memory
scampia_vis_wm = 3.5
roma_vis_wm = 3.5
pozzuoli_vis_wm = 3.8

#create the mean scores for each neighbourhood on each memory test 
hood_means = list(Scampia = c(verb_stm = scampia_verb_stm, 
                              vis_stm = scampia_vis_stm,
                              verb_wm = scampia_verb_wm,
                              vis_wm = scampia_vis_wm),
                  Roma = c(verb_stm = roma_verb_vstm,
                           vis_stm = roma_vis_vstm,
                           verb_wm = roma_verb_wm,
                           vis_wm = roma_vis_wm),
                  Pozzuoli = c(verb_stm = pozzuoli_verb_vstm,
                               vis_stm = pozzuoli_vis_vstm,
                               verb_wm = pozzuoli_verb_wm,
                               vis_wm = pozzuoli_vis_wm))

#create the standard deviation scores for each neighbourhood on each memory test 
hood_sds = list(Scampia = c(verb_stm = 0.75, 
                            vis_stm = 0.75,
                            verb_wm =  0.5,
                            vis_wm = 0.5),
                Roma = c(verb_stm = 0.75,
                         vis_stm = 0.75,
                         verb_wm = 0.5,
                         vis_wm = 0.5),
                Pozzuoli = c(verb_stm = 1,
                             vis_stm = 1,
                             verb_wm = 1,
                             vis_wm = 1))

#set the correlation of test score results for each neighbourhood
hood_cors = list(Scampia = .6, Roma = .6, Pozzuoli = .7)

#create the dataframe
dat = sim_design(within_short, between, n = 50, 
                 mu = hood_means, sd = hood_sds, r = hood_cors,
                 empirical = FALSE, plot = FALSE)

#make the data long format
dat_long = melt(setDT(dat), id.vars = c("id","neighbourhood"), variable.name = "test")

#ensure no values are below 2 or above 8
dat_long$value = ifelse(dat_long$value < 2, 2,
                        ifelse(dat_long$value > 8, 8, dat_long$value))

#set contrasts of neighbourhood so that "Pozzuoli" is the reference
contrasts(dat_long$neighbourhood) = contr.treatment(3, base = 3)

#round the test scores to a whole number
dat_long$value = round(dat_long$value, 0)

### ### ###

#plot the data and save it
test_labels = c("Verbal\nshort-term memory",
                "Visuospatial\nshort-term memory",
                "Verbal working\nmemory", 
                "Visuospatial\nworking memory")

scores_by_hood = ggplot(dat_long, aes(x = test, y = value)) + 
                  geom_violin(aes(fill = neighbourhood), trim = TRUE, position = position_dodge(0.9)) +
                  geom_boxplot(aes(fill = neighbourhood), width = 0.15, position = position_dodge(0.9)) +
                  scale_x_discrete(labels = test_labels) +
                  ylab("Test score") +
                  xlab("Test type") +
                  labs(fill = "Neighbourhood") +
                  aes(ymin = 0) +
                  avenir_theme

ggsave("../plots/plotted_simulated_data_study_1.jpg", scores_by_hood, width = 10, height = 5)

################################################################################################################################################

### ANALYSES ###

#create a dataset for each test type
verbstm_df = droplevels(subset(dat_long, dat_long$test == "verbstm"))
visstm_df = droplevels(subset(dat_long, dat_long$test == "visstm"))
verbwm_df = droplevels(subset(dat_long, dat_long$test == "verbwm"))
viswm_df = droplevels(subset(dat_long, dat_long$test == "viswm"))

### ### ###

#run models for each test type (first set the contrast for neighbourhood)
contrasts(verbstm_df$neighbourhood) = contr.treatment(3, base = 3)
verbstm_mod = lm(value ~ neighbourhood, data = verbstm_df)
summary(verbstm_mod)

### ### ###

library(effsize)

#create groups for each neighbourhood and test
roma_verbstm = subset(verbstm_df, verbstm_df$neighbourhood == "Roma")
roma_verbstm_list = as.numeric(roma_verbstm$value)

pozzuoli_verbstm = subset(verbstm_df, verbstm_df$neighbourhood == "Pozzuoli")
pozzuoli_verbstm_list = as.numeric(pozzuoli_verbstm$value)

#get Cohen's d for the comparisons between groups
cohen.d(roma_verbstm_list, pozzuoli_verbstm_list)

################################################################################################################################################

### DATA SIMULATION FOR POWER ANALYSIS ###

#create a function for power simulations (based on above data creation procedure); `participant_n` is the number of participants per neighbourhood
power_simulation = function(participant_n = 40,
                            pozzuoli_verb_vstm = 4.9,
                            pozzuoli_vis_vstm = 4.7,
                            pozzuoli_verb_wm = 3.9,
                            pozzuoli_vis_wm = 3.7,
                            high_stress_verb_effect = -0.4,
                            high_stress_vis_effect = -0.2,
                            pozzuoli_sd = 1,
                            high_stress_stm_sd = 0.75,
                            high_stress_wm_sd = 0.5,
                            test_score_cors_pozzuoli = 0.7,
                            test_score_cors_high_stress = 0.6,
                            ...) {
  
  #create the neighbourhood variable 
  between = list(neighbourhood = c(Scampia = "Scampia", 
                                   Roma = "Roma camp",
                                   Pozzuoli = "Pozzuoli "))
  
  #create the test variable
  within = list(test = c("Verbal short-term memory",
                         "Visuospatial short-term memory",
                         "Verbal working memory", 
                         "Visuospatial working memory"))
  
  within_short = list(test = c("verbstm",
                               "visstm",
                               "verbwm", 
                               "viswm"))
  
  #create scores for each neighbourhood for verbal short-term memory
  scampia_verb_stm = pozzuoli_verb_vstm + high_stress_verb_effect
  roma_verb_vstm = pozzuoli_verb_vstm + high_stress_verb_effect
  pozzuoli_verb_vstm = pozzuoli_verb_vstm
  
  #create scores for each neighbourhood for visuospatial short-term memory 
  scampia_vis_stm = pozzuoli_vis_vstm + high_stress_vis_effect
  roma_vis_vstm = pozzuoli_vis_vstm + high_stress_vis_effect
  pozzuoli_vis_vstm = pozzuoli_vis_vstm
  
  #create scores for each neighbourhood for verbal working memory
  scampia_verb_wm = pozzuoli_verb_wm + high_stress_verb_effect
  roma_verb_wm = pozzuoli_verb_wm + high_stress_verb_effect
  pozzuoli_verb_wm = pozzuoli_verb_wm
  
  #create scores for each neighbourhood for visuospatial working memory
  scampia_vis_wm = pozzuoli_vis_wm + high_stress_vis_effect
  roma_vis_wm = pozzuoli_vis_wm + high_stress_vis_effect
  pozzuoli_vis_wm = pozzuoli_vis_wm
  
  #create the mean scores for each neighbourhood on each memory test 
  hood_means = list(Scampia = c(verb_stm = scampia_verb_stm, 
                                vis_stm = scampia_vis_stm,
                                verb_wm = scampia_verb_wm,
                                vis_wm = scampia_vis_wm),
                    Roma = c(verb_stm = roma_verb_vstm,
                             vis_stm = roma_vis_vstm,
                             verb_wm = roma_verb_wm,
                             vis_wm = roma_vis_wm),
                    Pozzuoli = c(verb_stm = pozzuoli_verb_vstm,
                                 vis_stm = pozzuoli_vis_vstm,
                                 verb_wm = pozzuoli_verb_wm,
                                 vis_wm = pozzuoli_vis_wm))
  
  #create the standard deviation scores for each neighbourhood on each memory test 
  hood_sds = list(Scampia = c(verb_stm = high_stress_stm_sd, 
                              vis_stm = high_stress_stm_sd,
                              verb_wm =  high_stress_wm_sd,
                              vis_wm = high_stress_wm_sd),
                  Roma = c(verb_stm = high_stress_stm_sd,
                           vis_stm = high_stress_stm_sd,
                           verb_wm = high_stress_wm_sd,
                           vis_wm = high_stress_wm_sd),
                  Pozzuoli = c(verb_stm = pozzuoli_sd,
                               vis_stm = pozzuoli_sd,
                               verb_wm = pozzuoli_sd,
                               vis_wm = pozzuoli_sd))
  
  #set the correlation of test score results for each neighbourhood
  hood_cors = list(Scampia = test_score_cors_high_stress, Roma = test_score_cors_high_stress, Pozzuoli = test_score_cors_pozzuoli)
  
  #create the dataframe
  dat = sim_design(within_short, between, n = participant_n, 
                   mu = hood_means, sd = hood_sds, r = hood_cors,
                   empirical = FALSE, plot = FALSE)
  
  #make the data long format
  dat_long = melt(setDT(dat), id.vars = c("id","neighbourhood"), variable.name = "test")
  
  #ensure no values are below 2 or above 8
  dat_long$value = ifelse(dat_long$value < 2, 2,
                          ifelse(dat_long$value > 8, 8, dat_long$value))
  
  #set contrasts of neighbourhood so that "Pozzuoli" is the reference
  contrasts(dat_long$neighbourhood) = contr.treatment(3, base = 3)
  
  #round the test scores to a whole number
  dat_long$value = round(dat_long$value, 0)  
  
  ### ### ###

  #create a dataset for each test type
  verbstm_df = droplevels(subset(dat_long, dat_long$test == "verbstm"))
  visstm_df = droplevels(subset(dat_long, dat_long$test == "visstm"))
  verbwm_df = droplevels(subset(dat_long, dat_long$test == "verbwm"))
  viswm_df = droplevels(subset(dat_long, dat_long$test == "viswm"))
  
  ### ### ###
  
  #set the contrast for neighbourhood for each dataset
  contrasts(verbstm_df$neighbourhood) = contr.treatment(3, base = 3)
  contrasts(visstm_df$neighbourhood) = contr.treatment(3, base = 3)
  contrasts(verbwm_df$neighbourhood) = contr.treatment(3, base = 3)
  contrasts(viswm_df$neighbourhood) = contr.treatment(3, base = 3)
  
  #run model for each test type
  verbstm_mod = lm(value ~ neighbourhood, data = verbstm_df)
  visstm_mod = lm(value ~ neighbourhood, data = visstm_df)
  verbwm_mod = lm(value ~ neighbourhood, data = verbwm_df)
  viswm_mod = lm(value ~ neighbourhood, data = viswm_df)
  
  ### ### ###
  
  #create groups for one of the high-stress neighbourhoods (they are drawn from the same population) and calculate effect sizes for each outcome
  roma_verbstm = subset(verbstm_df, verbstm_df$neighbourhood == "Roma")
  roma_verbstm_list = as.numeric(roma_verbstm$value)
  pozzuoli_verbstm = subset(verbstm_df, verbstm_df$neighbourhood == "Pozzuoli")
  pozzuoli_verbstm_list = as.numeric(pozzuoli_verbstm$value)
  
  roma_visstm = subset(visstm_df, visstm_df$neighbourhood == "Roma")
  roma_visstm_list = as.numeric(roma_visstm$value)
  pozzuoli_visstm = subset(visstm_df, visstm_df$neighbourhood == "Pozzuoli")
  pozzuoli_visstm_list = as.numeric(pozzuoli_visstm$value)
  
  roma_verbwm = subset(verbwm_df, verbwm_df$neighbourhood == "Roma")
  roma_verbwm_list = as.numeric(roma_verbwm$value)
  pozzuoli_verbwm = subset(verbwm_df, verbwm_df$neighbourhood == "Pozzuoli")
  pozzuoli_verbwm_list = as.numeric(pozzuoli_verbwm$value)
  
  roma_viswm = subset(viswm_df, viswm_df$neighbourhood == "Roma")
  roma_viswm_list = as.numeric(roma_viswm$value)
  pozzuoli_viswm = subset(viswm_df, viswm_df$neighbourhood == "Pozzuoli")
  pozzuoli_viswm_list = as.numeric(pozzuoli_viswm$value)
  
  #get Cohen's d for the comparisons between groups
  effect_size_results_verbstm = cohen.d(roma_verbstm_list, pozzuoli_verbstm_list)
  effect_size_verbstm = as.numeric(effect_size_results_verbstm$estimate)
  effect_size_verbstm_lower = as.numeric(effect_size_results_verbstm$conf.int[1])
  effect_size_verbstm_upper = as.numeric(effect_size_results_verbstm$conf.int[2])
  
  effect_size_results_visstm = cohen.d(roma_visstm_list, pozzuoli_visstm_list)
  effect_size_visstm = as.numeric(effect_size_results_visstm$estimate)
  effect_size_visstm_lower = as.numeric(effect_size_results_visstm$conf.int[1])
  effect_size_visstm_upper = as.numeric(effect_size_results_visstm$conf.int[2])
  
  effect_size_results_verbwm = cohen.d(roma_verbwm_list, pozzuoli_verbwm_list)
  effect_size_verbwm = as.numeric(effect_size_results_verbwm$estimate)
  effect_size_verbwm_lower = as.numeric(effect_size_results_verbwm$conf.int[1])
  effect_size_verbwm_upper = as.numeric(effect_size_results_verbwm$conf.int[2])
  
  effect_size_results_viswm = cohen.d(roma_viswm_list, pozzuoli_viswm_list)
  effect_size_viswm = as.numeric(effect_size_results_viswm$estimate)
  effect_size_viswm_lower = as.numeric(effect_size_results_viswm$conf.int[1])
  effect_size_viswm_upper = as.numeric(effect_size_results_viswm$conf.int[2])
  
  ### ### ###

  #return a dataframe of the model results 
  verbstm_results = broom.mixed::tidy(verbstm_mod)
  verbstm_results$outcome = "Verbal short-term memory"
  verbstm_results$high_stress_effect = c(NA, effect_size_verbstm, NA)
  verbstm_results$high_stress_effect_lower = c(NA, effect_size_verbstm_lower, NA)
  verbstm_results$high_stress_effect_upper = c(NA, effect_size_verbstm_upper, NA)
  
  visstm_results = broom.mixed::tidy(visstm_mod)
  visstm_results$outcome = "Visuospatial short-term memory"
  visstm_results$high_stress_effect = c(NA, effect_size_visstm, NA)
  visstm_results$high_stress_effect_lower = c(NA, effect_size_visstm_lower, NA)
  visstm_results$high_stress_effect_upper = c(NA, effect_size_visstm_upper, NA)
  
  verbwm_results = broom.mixed::tidy(verbwm_mod)
  verbwm_results$outcome = "Verbal working memory"
  verbwm_results$high_stress_effect = c(NA, effect_size_verbwm, NA)
  verbwm_results$high_stress_effect_lower = c(NA, effect_size_verbwm_lower, NA)
  verbwm_results$high_stress_effect_upper = c(NA, effect_size_verbwm_upper, NA)
  
  viswm_results = broom.mixed::tidy(viswm_mod)
  viswm_results$outcome = "Visuospatial working memory"
  viswm_results$high_stress_effect = c(NA, effect_size_viswm, NA)
  viswm_results$high_stress_effect_lower = c(NA, effect_size_viswm_lower, NA)
  viswm_results$high_stress_effect_upper = c(NA, effect_size_viswm_upper, NA)
  
  do.call("rbind", list(verbstm_results, visstm_results, verbwm_results, viswm_results))

}

#run repeated simulations with dataset variants for each effect type
simulations_verbal = crossing(replications = 1:100,
                              participant_n = c(50, 75, 100, 125, 150),
                              pozzuoli_verb_vstm = 4.9,
                              pozzuoli_vis_vstm = 4.7,
                              pozzuoli_verb_wm = 3.9,
                              pozzuoli_vis_wm = 3.7,
                              high_stress_verb_effect = c(-0.3, -0.4, -0.5, -0.6, -0.7),
                              high_stress_vis_effect = -0.2,
                              pozzuoli_sd = 1,
                              high_stress_stm_sd = 0.75,
                              high_stress_wm_sd = 0.5,
                              test_score_cors_pozzuoli = 0.7,
                              test_score_cors_high_stress = 0.6) %>% mutate(analysis = pmap(., power_simulation)) %>% unnest(analysis)

simulations_visual = crossing(replications = 1:100,
                              participant_n = c(50, 75, 100, 125, 150),
                              pozzuoli_verb_vstm = 4.9,
                              pozzuoli_vis_vstm = 4.7,
                              pozzuoli_verb_wm = 3.9,
                              pozzuoli_vis_wm = 3.7,
                              high_stress_verb_effect = -0.4,
                              high_stress_vis_effect = c(-0.1, -0.2, -0.3, -0.4, -0.5),
                              pozzuoli_sd = 1,
                              high_stress_stm_sd = 0.75,
                              high_stress_wm_sd = 0.5,
                              test_score_cors_pozzuoli = 0.7,
                              test_score_cors_high_stress = 0.6) %>% mutate(analysis = pmap(., power_simulation)) %>% unnest(analysis)

#save the dataframes
write.csv(simulations_verbal, "../data/study1_verbal_memory_data_simulations.csv")
write.csv(simulations_visual, "../data/study1_visuospatial_memory_data_simulations.csv")

################################################################################################################################################

### POWER ANALYSIS ###

#subset the data to each outcome type
verb_stm_sims = subset(simulations_verbal, simulations_verbal$outcome == "Verbal short-term memory")
vis_stm_sims = subset(simulations_visual, simulations_visual$outcome == "Visuospatial short-term memory")
verb_wm_sims = subset(simulations_verbal, simulations_verbal$outcome == "Verbal working memory")
vis_wm_sims = subset(simulations_visual, simulations_visual$outcome == "Visuospatial working memory")

#create dataset to plot power analysis for neighbourhood main effects (both high-stress neighbourhoods were estimated to have the same effect)
simulation_results_neighbourhood_verb_stm = filter(verb_stm_sims, term == "neighbourhood1") %>%
                                            group_by(high_stress_verb_effect, participant_n) %>% 
                                            summarise(power = mean(p.value < .05), .groups = "drop")

simulation_results_neighbourhood_vis_stm = filter(vis_stm_sims, term == "neighbourhood1") %>%
                                           group_by(high_stress_vis_effect, participant_n) %>% 
                                           summarise(power = mean(p.value < .05), .groups = "drop")

simulation_results_neighbourhood_verb_wm = filter(verb_wm_sims, term == "neighbourhood1") %>%
                                           group_by(high_stress_verb_effect, participant_n) %>% 
                                           summarise(power = mean(p.value < .05), .groups = "drop")

simulation_results_neighbourhood_vis_wm = filter(vis_wm_sims, term == "neighbourhood1") %>%
                                           group_by(high_stress_vis_effect, participant_n) %>% 
                                           summarise(power = mean(p.value < .05), .groups = "drop")

### ### ###

#plot power analyse for neighbourhood effect on verbal short-term memory
verb_stm_pa = ggplot(aes(as.character(high_stress_verb_effect), participant_n, fill = power), 
                     data = simulation_results_neighbourhood_verb_stm) +
               geom_tile() +
               geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
               scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
               scale_fill_viridis_c(name = "Power",
                                    limits = c(0, 1), 
                                    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                    labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
               xlab("High-stress neighbourhood effect on verbal short-term memory") + 
               ylab("Participant sample size (per neighbourhood)") +
               avenir_theme

ggsave("../plots/power_analysis_verbal_short_term_memory.jpg", verb_stm_pa, width = 10, height = 5)

### ### ###

#plot power analyse for neighbourhood effect on visual short-term memory
vis_stm_pa = ggplot(aes(as.character(high_stress_vis_effect), participant_n, fill = power), 
                    data = simulation_results_neighbourhood_vis_stm) +
              geom_tile() +
              geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
              scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
              scale_fill_viridis_c(name = "Power",
                                   limits = c(0, 1), 
                                   breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                   labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
              xlab("High-stress neighbourhood effect on visual short-term memory") + 
              ylab("Participant sample size (per neighbourhood)") +
              avenir_theme

ggsave("../plots/power_analysis_visual_short_term_memory.jpg", vis_stm_pa, width = 10, height = 5)

### ### ###

#plot power analyse for neighbourhood effect on verbal working memory
verb_wm_pa = ggplot(aes(as.character(high_stress_verb_effect), participant_n, fill = power), 
                    data = simulation_results_neighbourhood_verb_wm) +
              geom_tile() +
              geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
              scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
              scale_fill_viridis_c(name = "Power",
                                   limits = c(0, 1), 
                                   breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                   labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
              xlab("High-stress neighbourhood effect on verbal working memory") + 
              ylab("Participant sample size (per neighbourhood)") +
              avenir_theme

ggsave("../plots/power_analysis_verbal_working_memory.jpg", verb_wm_pa, width = 10, height = 5)

### ### ###

#plot power analyse for neighbourhood effect on visual working memory
vis_wm_pa = ggplot(aes(as.character(high_stress_vis_effect), participant_n, fill = power), 
                   data = simulation_results_neighbourhood_vis_wm) +
              geom_tile() +
              geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
              scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
              scale_fill_viridis_c(name = "Power",
                                   limits = c(0, 1), 
                                   breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                   labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
              xlab("High-stress neighbourhood effect on visual working memory") + 
              ylab("Participant sample size (per neighbourhood)") +
              avenir_theme

ggsave("../plots/power_analysis_visual_working_memory.jpg", vis_wm_pa, width = 10, height = 5)

################################################################################################################################################

### EFFECT SIZE ESTIMATES ###

#create datasets
simulation_results_effect_sizes_verb_stm = verb_stm_sims[complete.cases(verb_stm_sims), ]
simulation_results_effect_sizes_vis_stm = vis_stm_sims[complete.cases(vis_stm_sims), ]

#select a participant sample size (per neighbourhood)
neighbourhood_n = 150
verb_stm_effect = -0.5
vis_stm_effect = -0.05

#get the mean effect size and CI for the assumed effect and a sample size of 150 for short-term memory (to estimate comparing effect sizes)
assumed_verb_stm_effect = subset(simulation_results_effect_sizes_verb_stm,
                                 simulation_results_effect_sizes_verb_stm$participant_n == neighbourhood_n & 
                                 simulation_results_effect_sizes_verb_stm$high_stress_verb_effect == verb_stm_effect)

assumed_vis_stm_effect = subset(simulation_results_effect_sizes_vis_stm,
                                simulation_results_effect_sizes_vis_stm$participant_n == neighbourhood_n & 
                                simulation_results_effect_sizes_vis_stm$high_stress_vis_effect == vis_stm_effect)


combined = data.frame(verb_stm_lower = assumed_verb_stm_effect$high_stress_effect_lower,
               verb_stm_effect = assumed_verb_stm_effect$high_stress_effect,
               verb_stm_upper = assumed_verb_stm_effect$high_stress_effect_upper,
               vis_stm_lower = assumed_vis_stm_effect$high_stress_effect_lower,
               vis_stm_effect = assumed_vis_stm_effect$high_stress_effect,
               vis_stm_upper = assumed_vis_stm_effect$high_stress_effect_upper)

#add a variable that is whether or not the confidence intervals overlap
combined$overlapping_CI = ifelse(combined$vis_stm_lower <= combined$verb_stm_upper, TRUE, FALSE)
print(paste0("PERCENTAGE OF SIGNIFICANT EFFECT SIZE DIFFERENCES: ",
            table(combined$overlapping_CI)[1] / sum(table(combined$overlapping_CI)) * 100, "%"))

print(paste0("MEAN EFEFCT SIZES: ", mean(assumed_verb_stm_effect$high_stress_effect)))
print(paste0("MEAN EFEFCT SIZES: ", mean(assumed_vis_stm_effect$high_stress_effect)))
