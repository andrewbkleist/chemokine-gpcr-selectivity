# Name:     7_tier3_quad_plot.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 4D

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ggrepel)
  library(bio3d)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")
  
##### ANALYSIS #################################################################
  
  # IMPORT DATA
  data <- read_csv("output/MASTER_CONSERVATION_WITH_TIERS.csv") 
  
  # IMPORT INTERFACE
  rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
    filter(Chain1 != Chain2) %>% select(target_gnccn) %>% distinct()
  rin$gn <- c("gn")
  inter <- unite(rin, gn, sep = "", c(gn, target_gnccn))
  inter <- inter$gn
  inter <- c(inter,  "1x28_", "1x24_", "1x25_") 
    # manually add positions adeed after sequence adjustment
  rm(rin)
 
  # MERGE
  data <- data %>% mutate(inter = case_when(
    gn %in% inter ~ "inter",
    !(gn %in% inter) ~ "no"
  ))
  
  # PLOT
  # ccl5 plotting
  ccr5 <- data %>% dplyr::filter(ckr == "ccr5") 
  ccr5 %>%
    ggplot(aes(cc, ortho_cons)) +
    geom_point(aes(fill = inter), shape = 21,size = 2, stroke = 0.5) +
    #geom_text_repel(data=subset(ccr5, cc < 0.8 & ortho_cons >= 0.8), aes(label=gn)) +
    scale_fill_manual(values=c('grey70', 'white'))+
    geom_hline(yintercept=0.8) +
    geom_vline(xintercept=0.8) +
    xlim(0,1) +
    theme_minimal()
  
  cxcr4 <- data %>% dplyr::filter(ckr == "cxcr4") 
  cxcr4 %>%
    ggplot(aes(cxc, ortho_cons)) +
    geom_point(aes(fill = inter), shape = 21,size = 2, stroke = 0.5) +
    #geom_text_repel(data=subset(ccr5, cc < 0.8 & ortho_cons >= 0.8), aes(label=gn)) +
    scale_fill_manual(values=c('grey70', 'white'))+
    geom_hline(yintercept=0.8) +
    geom_vline(xintercept=0.8) +
    xlim(0,1) +
    theme_minimal()
 