# Name:     6_ckr_core_logo.r
# Updated:  20191030
# User:     Andrew Kleist
# Paper:    Figure 4B

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  library(reshape2)
  library(RColorBrewer)
  require(Biostrings)
  library(ggseqlogo)
  library(seqRFLP)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")
  
##### 1: SEQUENCE LOGO #########################################################
  
  # load fasta files, convert to matrix, assign CCN
  receptor_master <- readAAMultipleAlignment("sequences/PARALOG_RECEPTOR_ALIGNMENT_CYS_ADJ.fasta")
  lookup_master <- as.matrix(receptor_master)
  lookup_master <- as.data.frame(lookup_master)
  
  gn_names <- c(read.table("sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt", sep = ",", colClasses = "character"))
  colnames(lookup_master) <- c(gn_names)
  
  
  # IMPORT TIER 1
  data <- read.csv("output/RECEPTOR_CONSERVATION.csv")
  tier1 <- data %>% filter(all_951 >= 0.8)
  tier1 <- as.character(unique(tier1$gn))
  tier1 <- c(tier1, "gn1x22_", "gn7x24") # manually add 1x22_ and 7x24
  
  # IMPORT TIER 2
  tier2 <- read_csv("output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(mean >= 0.8) # %>% filter(int == "int")
  tier2 <- tier2$motif
  tier1 <- c(tier1, "gn3x32") # manually add 3x32
  
  
  # IMPORT INTERFACE
  # import interface
  rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
    filter(Chain1 != Chain2) %>% select(target_gnccn) %>% distinct()
  rin$gn <- c("gn")
  inter <- unite(rin, gn, sep = "", c(gn, target_gnccn))
  inter <- inter$gn
  inter <- c(inter,  "gn1x28_", "gn1x24_", "gn1x25_")
  rm(rin)
  
  # SUBSET SEQUENCE POSITIONS BY TIER 1 AND TIER 2
  tier12 <- lookup_master[(names(lookup_master) %in% c(tier1, tier2))]
  
  # SUBSET SEQUENCE POSITIONS BY INTERFACE
  tier12int <- tier12[(names(tier12) %in% inter)]
  
  # MAKE SEQUENCE STRINGS
  core.all <- tier12[1:23,]
  core.all <- unite(core.all[,1:ncol(core.all)], col = seq_string,  sep = "")
  
  core.cc <- tier12[1:10,]
  core.cc <- unite(core.cc[,1:ncol(core.cc)], col = seq_string,  sep = "")
  
  core.cxc <- tier12[11:16,]
  core.cxc <- unite(core.cxc[,1:ncol(core.cxc)], col = seq_string,  sep = "")
  
  
  # GENERATE PLOTS
  p.all <- ggseqlogo(core.all,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p.cc <- ggseqlogo(core.cc,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p.cxc <- ggseqlogo(core.cxc,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(p.all, p.cc, p.cxc)
  #gridExtra::grid.arrange(p.cc, p.cxc)
  
  
  