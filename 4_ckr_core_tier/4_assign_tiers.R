# Name:     4_assign_tiers.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 4; Supp. Figure 8

# Note that this script performs tier assignment but does not plot the tier
# matrix in Supp. Figure 8

# Note that the table generated at the end has 1-7 positions per receptor with 
# NA tier values, corresponding to 9 distinct positions

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ggrepel)
  library(bio3d)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")

##### PART 1: SCORE TIERS 1, 2, & 3 ############################################
  
  # import data, etc
  data <- read.csv("output/RECEPTOR_CONSERVATION.csv")
  
  # add class information
  class <- read_csv("input/cc_cxc_lookup.csv")
  data$class <- class$cc_cxc[match(unlist(data$ckr), class$ckckr)]
  rm(class)
  
  # ASSIGN TIER 1 LABELS
  tier1 <- data %>% filter(all_951 >= 0.8)
  tier1 <- as.character(unique(tier1$gn))
  tier1 <- c(tier1, "gn1x22_", "gn7x24") # manually add 1x22_ and 7x24

  # IMPORT TIER 2 LABELS (get all tier 2 and subset later)
  tier2 <- read_csv("output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(mean >= 0.8) # %>% filter(int == "int")
  tier2 <- tier2$motif
    
  # CLEAN UP DATA FRAME
  # filter out ACK, etc
  data <- data %>% filter(class != "ack" & class != "xc" & class != "cx3c")
    
  # ADD TIERS
  data <- data %>% 
    dplyr::mutate(tier = dplyr::case_when(
      data$gn %in% tier1 & data$class == "cc" ~ "tier1",
      data$gn %in% tier1 & data$class == "cxc" ~ "tier1",
        
      data$gn %in% tier2 & data$class == "cc" ~ "tier2",
      data$gn %in% tier2 & data$class == "cxc" ~ "tier2",
        
      data$cc < 0.8 & data$ortho_cons >= 0.8 & data$class == "cc"  ~ "tier3",
      data$cxc < 0.8 & data$ortho_cons >= 0.8 & data$class == "cxc" ~ "tier3",
        
      data$ortho_cons < 0.8 & data$class == "cc" ~ "low_ortho",
      data$ortho_cons < 0.8 & data$class == "cxc" ~ "low_ortho"
        
      ))
    
  # override 3x32 as Tier2 (see methods)
  data <- data %>% mutate(tier = case_when(
    gn == "gn3x32" ~ "tier2",
    gn != "gn3x32" ~ tier
  ))
    
  # select core only
  core <- data %>% filter(dom != "NT" & dom != "CT" & dom != "ECL1" & dom != "ECL2" & dom != "ECL3" & dom != "ICL1" & dom != "ICL2" & dom != "ICL3")
  core2 <- c("gn23x50","gn23x51","gn23x52","gn34x50","gn34x51","gn34x52","gn34x53","gn34x54","gn34x55","gn34x56","gn34x57","gn34x58","gn45x50","gn45x51","gn45x52")
  core3 <- data %>% filter(gn %in% core2)
  core <- bind_rows(core, core3)
  rm(core2, core3)
    
  # overwrite "data"
  data <- core %>% distinct()
    
  # count total no. interface
  test <- data %>% dplyr::select(ckr, tier) %>% dplyr::count(ckr, tier)
  write_csv(data, "output/MASTER_CONSERVATION_WITH_TIERS.csv")
 
    x <- data %>% select(gn) %>% unique(
      
    )
    