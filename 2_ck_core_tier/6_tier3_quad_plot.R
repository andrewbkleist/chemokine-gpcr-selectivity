# Name:     6_tier3_quad_plot.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 2G 

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ape)
  library(stringr)
  library(ggrepel)
  library(seqinr)
  library(bio3d)
  library(seqRFLP)
  library(ggpubr)
  library(RColorBrewer)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/2_ck_core_tier/")

##### ANALYSIS ################################################
  
  # import data, as well as promiscuity data
  data <- read_csv("output/MASTER_CONSERVATION_WITH_TIERS.csv") 
  prom <- read_csv("input/no_ortho_20181129.csv")
  data <- left_join(data, prom)
  
  # add interface information
  int <- read_csv("input/ck_interface_lookup.csv")
  data <- left_join(data, int)
  data <- data %>% mutate(interface = case_when(
    interface == "int" ~ "int",
    is.na(interface) ~ "no"
  ))
  
  # ccl5 plotting
  ccl5 <- data %>% dplyr::filter(ck == "ccl5") 
  ccl5 %>%
    ggplot(aes(cc, ortho_cons)) +
    geom_point(aes(fill = interface), shape = 21,size = 2, stroke = 0.5) +
    geom_text_repel(data=subset(ccl5, cc<0.8 & ortho_cons>0.8), aes(label=ccn)) +
    scale_fill_manual(values=c('grey70', 'white')) +
    theme_minimal()
  
