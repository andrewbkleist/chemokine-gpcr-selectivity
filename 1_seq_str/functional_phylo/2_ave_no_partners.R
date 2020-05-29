# Name:     2_ave_no_partners.R
# Updated:  20190401
# User:     Andrew Kleist
# Figure:   Figure 1G

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
  ck <- ck.fnc %>% count(ck)
  ckr <- ck.fnc %>% count(ckr)
  colnames(ck)[1] <- c("protein")
  colnames(ckr)[1] <- c("protein")
  ck$type <- c("ck")
  ckr$type <- c("ckr")
  
  # add in chemokines for which number of contacts equal zero
  # including CCL4L1, CXCL17, CXLC14, and CXCL4L1
  missing <-  data.frame(protein = c("CCL4L1", "CXCL4L1", "CXCL17", "CXCL14"), 
                         n = c(0,0,0,0),  type = c("ck", "ck", "ck", "ck"))
  ck <- bind_rows(ck, missing)
  rm(missing)
  
  # VIOLON PLOT - CHEMOKINES VS. RECEPTORS
  master <- bind_rows(ck, ckr)
  master %>% 
    ggplot(aes(type, n)) +
    geom_violin(trim=TRUE) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1, fill = "white", color="grey50") +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
    
    #geom_boxplot(width=0.05, bins=5, outlier.shape = NA) +
    theme_minimal()
  
  # GET MEANS AND MEDIANS
  mean(ck$n)
  mean(ckr$n)
  median(ck$n)
  median(ckr$n)
  
  # VIOLON PLOT - CHEMOKINES VS. RECEPTORS CC vs CXC
  lookup <- read_csv("input/cc_cxc_lookup.csv")
  lookup$ckckr <- toupper(lookup$ckckr)
  master$cc_cxc <- lookup$cc_cxc[match(unlist(master$protein), lookup$ckckr)]

  master %>% 
    filter(cc_cxc != "ack") %>%
    filter(cc_cxc != "xc") %>%
    filter(cc_cxc != "cx3c") %>%
    ggplot(aes(cc_cxc, n)) +
    geom_violin(trim=FALSE) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1, fill = "white", color="grey50") +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
    facet_grid(. ~ type) +
    
    #geom_boxplot(width=0.05, bins=5, outlier.shape = NA) +
    theme_minimal()
    
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
  
  
  
  
  