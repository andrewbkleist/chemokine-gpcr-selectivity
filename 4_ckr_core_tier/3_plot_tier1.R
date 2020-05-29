# Name:     3_plot_tier1.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 7 (left panel)

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")

##### 1: TIER 1 PLOTTING #######################################################
  
  # import
  data <- read_csv("output/RECEPTOR_CONSERVATION.csv")
  
  # remove gn name
  data <- data %>% separate(gn, c("gn1", "gn"), sep = "gn") %>% select(-gn1)
  
  # select all sequences
  data <- data %>% select(gn, all_951) %>% distinct()
  
  # import interface
  rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
    filter(Chain1 != Chain2) %>% select(target_gnccn) %>% distinct()
  inter <- rin$target_gnccn
  rm(rin)
  
  # order and plot
  data$gn <- factor(data$gn, levels = data$gn[order(data$all_951)])
  data %>% #filter(all_951 >= 0.8) %>%
    top_n(50, all_951) %>%
    #filter(gn %in% inter) %>%
    #filter(gn != "gn45x50" & gn != "gn45x51" & gn != "gn45x52") %>%
    ggplot(aes(gn, all_951)) +
    geom_point() +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()  
  
