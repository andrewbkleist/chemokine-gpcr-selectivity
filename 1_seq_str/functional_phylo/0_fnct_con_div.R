# Name:     0_fnct_con_div.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 1E,F 

##### LOAD PACKAGES & SET WD ###################################################

  # load libraries, set wd
  library(tidyverse)
  library(ape)
  library(phangorn)
  require(Biostrings)
  library(ggrepel)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/1_seq_str/functional_phylo/")
  
  
##### PART 1: CHEMOKINE - FUNCTIONAL DIVERGENCE CONVERGENECE SCATTER ###########
  
  # (1.1) MAKE FUNCTIONAL DISTANCE MATRIX
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
    
    # calculate distance
    ck.dist <- dist(ck.fnc, method = "euclidean")
    ck.dist <- as.data.frame(as.matrix(ck.dist))
    ck.dist$ck1 <- rownames(ck.dist)
    ck.dist <- ck.dist %>% gather(ck2, fnc_dist, 1:46)
    ck.dist <- ck.dist %>% mutate(fnc_dist_nl = fnc_dist / max(fnc_dist))
    
    # REMOVE UNDERSTUDIED CHEMOKINES
    under <- c("CCL3L1","CCL4L1","CXCL4L1","CXCL17","CXCL14") # remove chemokines with no reported interactions over 7
    ck.dist <- ck.dist %>% filter(!ck1 %in% under)
    ck.dist <- ck.dist %>% filter(!ck2 %in% under)
    
    # unite name
    ck.dist <- unite_(ck.dist, "ck1ck2", c("ck1", "ck2"))
    
    # remove used objects
    rm(ck.fnc)
    
  # (1.2) SEQUENCE DISTANCE MATRIX
    # Make sequence x sequence identity
    # Ran unaligned core sequences through ClustalO:
    # https://www.ebi.ac.uk/Tools/msa/clustalo/
    # took output pairwise seq identities
    ck.seq <- read.csv("input/identity_matrix_mod.csv", header = FALSE)
    names <- ck.seq$V1
    rownames(ck.seq) <- ck.seq[,1]
    ck.seq <- ck.seq[,-1]
    colnames(ck.seq) <- rownames(ck.seq)
    
    # gather,format
    ck.seq$ck1 <- names
    ck.seq <- ck.seq %>% gather(key = "ck2", value = "seq_id", 1:46)
    ck.seq <- unite_(ck.seq, "ck1ck2", c("ck1", "ck2"))
    ck.seq <- ck.seq %>% mutate(seq_id = seq_id /100)
    ck.seq <- ck.seq %>% mutate(seq_dist = 1- seq_id)
    
    # REMOVE UNDERSTUDIED CHEMOKINES
    under <- c("CCL3L1","CCL4L1","CXCL4L1","CXCL17","CXCL14") # remove chemokines with no reported interactions over 7
    ck.dist <- ck.dist %>% filter(!ck1 %in% under)
    ck.dist <- ck.dist %>% filter(!ck2 %in% under)
    
    # combine
    seq.fnc <- left_join(ck.dist, ck.seq)
    
    # remove used
    rm(ck.dist, ck.seq)
    
  # (1.3) PLOT
    seq.fnc %>%
      filter(seq_id !=1) %>%
      ggplot(aes(seq_id, fnc_dist_nl)) +
      geom_point(shape = 21, colour = "grey30", fill = "white", size = 2, stroke = 0.5) +
      xlim(0,1)+
      ylim(0,1)+
      #geom_text_repel(data=subset(seq.fnc, seq_dist>0.85 & fnc_dist_nl>0.85), aes(label=ck1ck2)) +
      theme_minimal()
    seq.fnc %>%
      filter(ck1ck2 == "CXCL9_CXCL10") %>%
      ggplot(aes(seq_id, fnc_dist_nl)) +
      geom_point(shape = 21, colour = "grey30", fill = "white", size = 2, stroke = 0.5) +
      xlim(0,1)+
      ylim(0,1)+
      #geom_text_repel(data=subset(seq.fnc, seq_dist>0.6 & fnc_dist_nl<0.3), aes(label=ck1ck2)) +
      theme_minimal()
  
    
  # (1.4) PIE CHART
    seq.fnc <- seq.fnc %>% filter(seq_id !=1)
    seq.fnc <- seq.fnc %>% mutate(quad = case_when(
      fnc_dist_nl > 0.5 & seq_id < 0.5 ~ "TR",
      fnc_dist_nl > 0.5 & seq_id >= 0.5 ~ "div",
      fnc_dist_nl < 0.5 & seq_id < 0.5 ~ "conv",
      fnc_dist_nl < 0.5 & seq_id >= 0.5 ~ "simsim"
    ))
    fnc.sum <- seq.fnc %>% count(quad)
    fnc.sum <- fnc.sum %>% mutate(pct = n / sum(n)*100)
      
    fnc.sum %>%
      ggplot(aes(x = "", y= pct, fill = quad)) +
      geom_bar(width = 1,size = 1, stat="identity", color = "white") +
      coord_polar("y") +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())


##### PART 2: RECEPTOR - FUNCTIONAL DIVERGENCE CONVERGENECE SCATTER ############
    
    # (1.1) MAKE FUNCTIONAL DISTANCE MATRIX
    fnc <- read.csv("input/supp_table_1.csv") # import functional table
    fnc[is.na(fnc)] <- 0 # unstudied interactions with zero
    
    fnc <- fnc %>% select(-H4.Receptor, -CXCR8...GPR35)
    rownames(fnc) <- fnc[,1] # replace row names w chemokines
    fnc <- fnc[,-1] # remove chemokines as row
      # this analysis does not take into account agonist vs. antagonist
    fnc <- t(fnc)
    fnc <- as.data.frame(fnc)
    fnc <- fnc %>% select(-CXCL15)
    
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
    
    # calculate distance
    dist <- dist(fnc, method = "euclidean")
    dist <- as.data.frame(as.matrix(dist))
    dist$ckr1 <- rownames(dist)
    dist <- dist %>% gather(ckr2, fnc_dist, 1:23)
    dist <- dist %>% mutate(fnc_dist_nl = fnc_dist / max(fnc_dist))

    # unite name
    dist <- unite_(dist, "ckr1ckr2", c("ckr1", "ckr2"))
    
    # remove used objects
    rm(fnc)
    
    # SEQUENCE
    # Make sequence x sequence identity
    # Ran unaligned core sequences through ClustalO:
    # https://www.ebi.ac.uk/Tools/msa/clustalo/
    # took output pairwise seq identities
    all.ckr.seq <- read.table("input/ckr_clustal_pid_matrix.csv", header = FALSE)
    names <- all.ckr.seq[,1]
    rownames(all.ckr.seq) <- all.ckr.seq[,1]
    all.ckr.seq <- all.ckr.seq[,-1]
    colnames(all.ckr.seq) <- rownames(all.ckr.seq)
    
    # gather, unite
    all.ckr.seq$ckr1 <- names
    all.ckr.seq <- gather(all.ckr.seq, key = "ckr2", value = "pid", 1:23)
    all.ckr.seq <- unite_(all.ckr.seq, "ckr1ckr2", c("ckr1", "ckr2"))
    all.ckr.seq$dist <- 100 - all.ckr.seq$pid  

    # combine
    seq.fnc <- left_join(dist, all.ckr.seq)
    
    # make sequence distance 0-to-1 scale
    seq.fnc <- seq.fnc %>% mutate(dist = dist/100)
    
    # remove used
    rm(dist, all.ckr.seq)
    
    # (1.3) PLOT
    seq.fnc %>%
      filter(pid !=100) %>%
      ggplot(aes(pid/100, fnc_dist_nl)) +
      geom_point(shape = 21, colour = "grey30", fill = "white", size = 2, stroke = 0.5) +
      xlim(0,1)+
      ylim(0,1)+
      #geom_text_repel(data=subset(seq.fnc, seq_dist>0.85 & fnc_dist_nl>0.85), aes(label=ck1ck2)) +
      theme_minimal()
    seq.fnc %>%
      filter(ckr1ckr2 == "CXCR1_CXCR2") %>%
      ggplot(aes(pid/100, fnc_dist_nl)) +
      geom_point(shape = 21, colour = "grey30", fill = "white", size = 2, stroke = 0.5) +
      xlim(0,1)+
      ylim(0,1)+
      #geom_text_repel(data=subset(seq.fnc, seq_dist>0.6 & fnc_dist_nl<0.3), aes(label=ck1ck2)) +
      theme_minimal()
    
    # (1.5) PIE CHART
    seq.fnc <- seq.fnc %>% filter(dist !=0)
    
    seq.fnc <- seq.fnc %>% mutate(quad = case_when(
      fnc_dist_nl > 0.5 & pid > 50 ~ "TR",
      fnc_dist_nl > 0.5 & pid <= 50 ~ "div",
      fnc_dist_nl < 0.5 & pid > 50 ~ "conv",
      fnc_dist_nl < 0.5 & pid <= 50 ~ "simsim"
    ))
    fnc.sum <- seq.fnc %>% count(quad)
    fnc.sum <- fnc.sum %>% mutate(pct = n / sum(n)*100)
    
    fnc.sum %>%
      ggplot(aes(x = "", y= pct, fill = quad)) +
      geom_bar(width = 1,size = 1, stat="identity", color = "white") +
      coord_polar("y") +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  