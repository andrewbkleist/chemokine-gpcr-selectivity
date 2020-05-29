# Name:     6_functional_valid.R
# Updated:  20191031
# User:     Andrew Kleist
# Figure:   Figure 3I

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/20190201_FINAL/3_ck_motif/")

##### GET FUNCTIONAL ##############################################################

  # cxcl8
  data <- read_csv("input/elr_mutations.csv")
  
  order <- c("wt", "4_72", "5_72", "6_72", "7_72")
  data$mutant <- factor(data$mutant, levels = rev(order))
  
  data %>%
    ggplot(aes(mutant, kdrat)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    ylim(0,220) +
    coord_flip()
  
  # ccl5
  data <- read_csv("input/ccl5_mutations.csv")
  data <- data %>% mutate(ec50rat = ec50/29.1)
  
  order <- c("wt", "3_68", "4_68")
  data$mutant <- factor(data$mutant, levels = rev(order))
  
  data %>%
    ggplot(aes(mutant, ec50rat)) +
    geom_bar(stat="identity") +
    # geom_errorbar(aes(ymin=ec50-error, ymax=ec50+error), width=.2,
    #               position=position_dodge(.9)) +
    theme_minimal() +
    #ylim(0,220) +
    coord_flip()
  