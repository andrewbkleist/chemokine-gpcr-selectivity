# Name:     5_ckr_core_logo.r
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 2F


##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  library(reshape2)
  library(RColorBrewer)
  require(Biostrings)
  library(ggseqlogo)
  library(seqRFLP)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/2_ck_core_tier/")

##### ANALYSIS #################################################################
  
  # load fasta files, convert to matrix, assign CCN
  master <- readAAMultipleAlignment("sequences/ALL_para.fasta")
  lookup_master <- as.matrix(master)
  lookup_master <- as.data.frame(lookup_master)
  
  gn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
  colnames(lookup_master) <- c(gn_names)
  
  
  # IMPORT TIER 1
  data <- read.csv("output/CHEMOKINE_CONSERVATION.csv")
  tier1 <- c("CX.1", "CX.5", "b1b2.12", "B3.3", "H.3", "b3h.2")
  
  # IMPORT TIER 2
  tier2 <- read_csv("output/CHEMOKINE_CLASSIFICATION_N3.csv") %>%
    filter(mean >= 0.8) # %>% filter(int == "int")
  tier2 <- tier2$motif
  
  # IMPORT INTERFACE
  # import interface
  rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
    filter(Chain1 != Chain2) %>% select(source_gnccn) %>% distinct()
  inter <- rin$source_gnccn
  rm(rin)
  
  # SUBSET SEQUENCE POSITIONS BY TIER 1 AND TIER 2
  tier12 <- lookup_master[(names(lookup_master) %in% c(tier1, tier2))]
  
  # SUBSET SEQUENCE POSITIONS BY INTERFACE
  tier12int <- tier12[(names(tier12) %in% inter)]
  
  # MAKE SEQUENCE STRINGS
  core.all <- tier12int[1:43,]
  core.all <- unite(core.all[,1:ncol(core.all)], col = seq_string,  sep = "")
  
  core.cc <- tier12[1:26,]
  core.cc <- unite(core.cc[,1:ncol(core.cc)], col = seq_string,  sep = "")
  
  core.cxc <- tier12[27:43,]
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
  
  
  