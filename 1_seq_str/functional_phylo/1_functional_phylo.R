# Name:     1_seq_struct_phylo.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 1A,D


##### LOAD PACKAGES & SET WD ###################################################

  # load libraries and wd
  library(tidyverse)
  library(ape)
  library(phangorn)
  require(Biostrings)
  library(ggrepel)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/1_seq_str/functional_phylo/")
  
  
##### NETWORK AND CLUSTERING TREE ##############################################
  
  # (1.1) FUNCTIONAL CLUSTERING TREE - CHEMOKINE -------------------------------
  
    ck.fnc <- read.csv("input/supp_table_1.csv") # import functional table
    ck.fnc[is.na(ck.fnc)] <- 0 # unstudied interactions with zero
    ck.fnc <- ck.fnc %>% filter(X != "CXCL15")
    rownames(ck.fnc) <- ck.fnc[,1] # replace row names w chemokines
    ck.fnc <- ck.fnc[,-1] # remove chemokines as row
    ck.fnc <- ck.fnc %>% dplyr::select(-H4.Receptor, -CXCR8...GPR35) # remove other receptors
    # this analysis does not take into account agonist vs. antagonist
    
    # change to 0 to 5 scale
    ck.fnc$ck <- rownames(ck.fnc)
    ck.fnc <- ck.fnc %>% gather(ckr, score, 1:23)
    ck.fnc <- ck.fnc %>% mutate(score = case_when(
      score == 10 ~ 5,
      score == 9 ~ 5,
      score == 8 ~ 4,
      score == 7 ~ 4,
      score == 6 ~ 3,
      score == 5 ~ 3,
      score == 4 ~ 2,
      score == 3 ~ 2,
      score == 2 ~ 1,
      score == 1 ~ 1,
      score == 0 ~ 0
    ))
    ck.fnc <- ck.fnc %>% spread(ckr, score)
    rownames(ck.fnc) <- ck.fnc$ck
    ck.fnc <- ck.fnc %>% select(-ck)
    
    # calculate distances
    ck.dist <- dist(ck.fnc, method = "euclidean")
    hc <- hclust(ck.dist, method = "ward.D2")  
    plot(hc,  hang = -1, cex = 0.6)
    # plot(as.phylo(hc),  hang = -1, cex = 0.6, type = "fan")
    
    # remove used objects
    rm(hc, ck.dist)
    
  # (1.2) FUNCTIONAL CLUSTERING TREE - RECEPTOR --------------------------------
    
    ck.fnc <- read.csv("input/supp_table_1.csv") # import functional table
    ck.fnc[is.na(ck.fnc)] <- 0 # unstudied interactions with zero
    ck.fnc <- ck.fnc %>% filter(X != "CXCL15")
    rownames(ck.fnc) <- ck.fnc[,1] # replace row names w chemokines
    ck.fnc <- ck.fnc[,-1] # remove chemokines as row
    ck.fnc <- ck.fnc %>% dplyr::select(-H4.Receptor, -CXCR8...GPR35) # remove other receptors
    # this analysis does not take into account agonist vs. antagonist
    
    # change to 0 to 5 scale
    ck.fnc$ck <- rownames(ck.fnc)
    ck.fnc <- ck.fnc %>% gather(ckr, score, 1:23)
    ck.fnc <- ck.fnc %>% mutate(score = case_when(
      score == 10 ~ 5,
      score == 9 ~ 5,
      score == 8 ~ 4,
      score == 7 ~ 4,
      score == 6 ~ 3,
      score == 5 ~ 3,
      score == 4 ~ 2,
      score == 3 ~ 2,
      score == 2 ~ 1,
      score == 1 ~ 1,
      score == 0 ~ 0
    ))
    ck.fnc <- ck.fnc %>% spread(ckr, score)
    rownames(ck.fnc) <- ck.fnc$ck
    ck.fnc <- ck.fnc %>% select(-ck)
    
    # transpose
    ckr.fnc <- t(ck.fnc)
    
    ckr.dist <- dist(ckr.fnc, method = "euclidean")
    hc <- hclust(ckr.dist, method = "ward.D2")  
    plot(hc,  hang = -1, cex = 0.6)
    # plot(as.phylo(hc),  hang = -1, cex = 0.6, type = "fan")
    
    # remove used objects
    rm(hc, ck.dist)
    
  # (1.2) CHEMOKINE-GPCR COUPLING MATRIX ---------------------------------------
    
    # gather for heatmap
    ck.fnc$ck <- rownames(ck.fnc)
    ck.fnc <- ck.fnc %>% gather(ckr, score, 1:(ncol(ck.fnc)-1))
    
    # reorder for heatmap
    order.ck <- toupper(c("ccl1","ccl2","ccl3","ccl3l1","ccl4","ccl4l1","ccl5","ccl6","ccl7","ccl8","ccl9","ccl10","ccl11","ccl12","ccl13","ccl14","ccl15","ccl16","ccl17","ccl18","ccl19","ccl20","ccl21","ccl22","ccl23","ccl24","ccl25","ccl26","ccl27","ccl28","cxcl1","cxcl2","cxcl3","cxcl4","cxcl4l1","cxcl5","cxcl6","cxcl7","cxcl8","cxcl9","cxcl10","cxcl11","cxcl12","cxcl13","cxcl14","cxcl15","cxcl16","cxcl17","cx3cl1","xcl1","xcl2")) 
    order.ckr <- toupper(c("ccr1","ccr2","ccr3","ccr4","ccr5","ccr6","ccr7","ccr8","ccr9","ccr10","cxcr1","cxcr2","cxcr3","cxcr4","cxcr5","cxcr6", "cxcr8", "xcr1","cx3cr1","cx3c1","ackr1","ackr2","ackr3","ackr4","ccrl2"))
    
    ck.fnc$ck <- factor(ck.fnc$ck, levels = rev(order.ck))
    ck.fnc$ckr <- factor(ck.fnc$ckr, levels = order.ckr)
    
    # plot
    ck.fnc %>%
      ggplot() +
      geom_tile(aes(ckr, ck, fill = score)) +
      #coord_fixed() +
      theme_minimal() +
      scale_fill_gradient(low="white", high="grey40") +
      #coord_flip() +
      theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

 
    