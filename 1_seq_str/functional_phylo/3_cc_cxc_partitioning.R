# Name:     3_ave_no_partners.R
# Updated:  20200406
# User:     Andrew Kleist

##### LOAD PACKAGES & SET WD ###################################################

  # load libraries and wd
  library(tidyverse)
  library(ape)
  library(phangorn)
  require(Biostrings)
  library(ggrepel)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/1_seq_str/functional_phylo/")
  
  
##### ANALYSIS #################################################################
  
  # import, reformat
  ck.fnc <- read.csv("input/supp_table_1.csv") # import functional table
  ck.fnc[is.na(ck.fnc)] <- 0 # unstudied interactions with zero
  ck.fnc <- ck.fnc %>% filter(X != "CXCL15")
  rownames(ck.fnc) <- ck.fnc[,1] # replace row names w chemokines
  colnames(ck.fnc)[1] <- c("ck")
  rownames(ck.fnc) <- 1:nrow(ck.fnc)
  ck.fnc <- ck.fnc %>% dplyr::select(-H4.Receptor, -CXCR8...GPR35) # remove other receptors
  ck.fnc <- ck.fnc %>% gather(ckr, score, CCR1:CCRL2)
  
  # subset, count - ONLY INTERACTIONS WITH SCORES > 6 ARE COUNTED
  ck.fnc <- ck.fnc %>% filter(score > 6) %>% select(-score)
  
  # add cc cxc
  lookup <- read_csv("input/cc_cxc_lookup.csv")
  lookup$ckckr <- toupper(lookup$ckckr)
  ck.fnc$cc_cxc_c <- lookup$cc_cxc[match(unlist(ck.fnc$ck), lookup$ckckr)]
  ck.fnc$cc_cxc_r <- lookup$cc_cxc[match(unlist(ck.fnc$ckr), lookup$ckckr)]
  
  # count number of cc-to-cc, cc-to-cxc, etc
  ck.fnc <- ck.fnc %>% count(cc_cxc_c, cc_cxc_r)
  ck.fnc <- ck.fnc %>% filter(cc_cxc_c == "cc" | cc_cxc_c == "cxc" | cc_cxc_c == "ack") %>%
    filter(cc_cxc_r == "cc" | cc_cxc_r == "cxc" | cc_cxc_r == "ack")
  
  test <- data.frame("CC" = c(50,0) , "CXC" = c(5,18))
  colnames(test) <- c("CC", "CXC")
  chisq.test(test)
  
  
  
  # GET MEANS AND MEDIANS
  master %>% 
    filter(cc_cxc == "cc") %>%
    filter(type == "ck") %>%
    summarise(mean = mean(n))
  
  master %>% 
    filter(cc_cxc == "cxc") %>%
    filter(type == "ck") %>%
    summarise(mean = mean(n))
  
  master %>% 
    filter(cc_cxc == "cc") %>%
    filter(type == "ckr") %>%
    summarise(mean = median(n))
  
  master %>% 
    filter(cc_cxc == "cxc") %>%
    filter(type == "ckr") %>%
    summarise(mean = median(n))
  
  # STATS
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  ck.test <-  master %>%
    filter(cc_cxc != "ack") %>%
    filter(cc_cxc != "cx3c") %>%
    filter(cc_cxc != "xc") %>%
    filter(type == "ck") 
  
  ck.test <- wilcox.test(n ~ cc_cxc, data = ck.test, exact = FALSE)
  ck.test$p.value
 
  ckr.test <-  master %>%
    filter(cc_cxc != "ack") %>%
    filter(cc_cxc != "cx3c") %>%
    filter(cc_cxc != "xc") %>%
    filter(type == "ckr") 
  
  ckr.test <- wilcox.test(n ~ cc_cxc, data = ckr.test, exact = FALSE)
  ckr.test$p.value
  
  # NETWORK PARTITIONING
  
  
  
  
  