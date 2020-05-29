# Name:     5_arg_vs_asp.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 5D

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/5_ckr_motif/")

##### GET CHARGED NTERM ##############################################################

  # (1) PREP DATA FRAME
  # import, reformat
  data <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")  

  # select only CONSERVED MOTIFS
  data <- data %>% filter(pct_ortho >= 0.5)

  # (3) ADD LABELS
  cc.cxc.ack <- read.csv("input/cc_cxc_ack.csv")
  data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
  rm(cc.cxc.ack)
  
  # (4) COUNT HOW MANY CC AND CXC MOTIFS you are querying from FOR NORMALIZATION
  how.many <- data %>% filter(pct_ortho >= 0.5)
  how.many <- how.many %>% select(class) %>% count(class) %>% 
    filter(class == "cc" | class == "cxc")
  print(how.many) # 621 CC, 356 CXC
  
  # (5) COUNT HOW MANY MOTIFS CONTAIN CHARGED RESIDUES
  asp <- subset(data, grepl("D", data$motif))
  asp$aa <- c("asp")
  glu <- subset(data, grepl("E", data$motif))
  glu$aa <- c("glu")
  arg <- subset(data, grepl("R", data$motif))
  arg$aa <- c("arg")
  lys <- subset(data, grepl("K", data$motif))
  lys$aa <- c("lys")
  his <- subset(data, grepl("H", data$motif))
  his$aa <- c("his")
  
  # (6) COUNT HOW MANY MOTIFS CONTAIN TYROSINES
  tyr <- subset(data, grepl("Y", data$motif))
  tyr$aa <- c("tyr")

  # (7) BIND CHARGED MOTIFS TOGETHER & ADD PLUS MINUS
  charged <- bind_rows(asp, glu, lys, arg, his)
  #rm(asp, glu, arg, lys, his)
  
  # add plus minus
  charged <- charged %>% mutate(charge = case_when(
    aa == "asp" ~ "minus",
    aa == "glu" ~ "minus",
    aa == "arg" ~ "plus",
    aa == "lys" ~ "plus",
    aa == "his" ~ "plus"
  ))
  
  # remove residues counted twice
  charged <- charged %>% select(motif, protein, class, charge)
  charged <- charged %>% distinct()
  
  # (8) NORMALIZE BY NUMBER OF RESIDUES QUERIED WITHIN EACH FAMILY
  class <- charged %>% count(class)
  charged <- charged %>%
    count(charge, class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/621,
      class == "cxc" ~ n/356
    ))
  
  # (9) BAR PLOT - CHARGE
  charged %>%
    filter(class != "cx3c" & class != "xc" & class != "ack") %>%
    ggplot(aes(charge, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.6) +
    theme_minimal() +
    facet_grid(. ~ class)
  
  # (10) TYROSINE NORMALIZATION
  class <- tyr %>% count(class)
  tyr <- tyr %>%
    count(class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/621,
      class == "cxc" ~ n/356
    ))

  # (11) BAR PLOT - TYROSINE
  tyr %>%
    filter(class != "cx3c" & class != "xc" & class != "ack") %>%
    ggplot(aes(class, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.4) +
    theme_minimal()
  
##### GET CHARGED ECL2 ##############################################################
  
  # (1) PREP DATA FRAME
  # import, reformat
  data <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv")  
  
  # select only CONSERVED MOTIFS
  data <- data %>% filter(pct_ortho >= 0.5)
  
  # (3) ADD LABELS
  cc.cxc.ack <- read.csv("input/cc_cxc_ack.csv")
  data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
  rm(cc.cxc.ack)
  
  # (4) COUNT HOW MANY CC AND CXC MOTIFS you are querying from FOR NORMALIZATION
  how.many <- data %>% filter(pct_ortho >= 0.5)
  how.many <- how.many %>% select(class) %>% count(class) %>% 
    filter(class == "cc" | class == "cxc")
  print(how.many) # 288 CC, 169 CXC
  
  # (5) COUNT HOW MANY MOTIFS CONTAIN CHARGED RESIDUES
  asp <- subset(data, grepl("D", data$motif))
  asp$aa <- c("asp")
  glu <- subset(data, grepl("E", data$motif))
  glu$aa <- c("glu")
  arg <- subset(data, grepl("R", data$motif))
  arg$aa <- c("arg")
  lys <- subset(data, grepl("K", data$motif))
  lys$aa <- c("lys")
  his <- subset(data, grepl("H", data$motif))
  his$aa <- c("his")
  
  # (6) COUNT HOW MANY MOTIFS CONTAIN TYROSINES
  tyr <- subset(data, grepl("Y", data$motif))
  tyr$aa <- c("tyr")
  
  # (7) BIND CHARGED MOTIFS TOGETHER & ADD PLUS MINUS
  charged <- bind_rows(asp, glu, lys, arg, his)
  #rm(asp, glu, arg, lys, his)
  
  # add plus minus
  charged <- charged %>% mutate(charge = case_when(
    aa == "asp" ~ "minus",
    aa == "glu" ~ "minus",
    aa == "arg" ~ "plus",
    aa == "lys" ~ "plus",
    aa == "his" ~ "plus"
  ))
  
  # remove residues counted twice
  charged <- charged %>% select(motif, protein, class, charge)
  charged <- charged %>% distinct()
  
  # (8) NORMALIZE BY NUMBER OF RESIDUES QUERIED WITHIN EACH FAMILY
  class <- charged %>% count(class)
  charged <- charged %>%
    count(charge, class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/288,
      class == "cxc" ~ n/169
    ))
  
  # (9) BAR PLOT - CHARGE
  charged %>%
    filter(class != "cx3c" & class != "xc" & class != "ack") %>%
    ggplot(aes(charge, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.6) +
    theme_minimal() +
    facet_grid(. ~ class)
  
  
##### GET AROMATIC ECL2 ########################################################
  
  # (1) PREP DATA FRAME
  # import, reformat
  data <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv")  
  
  # select only CONSERVED MOTIFS
  data <- data %>% filter(pct_ortho >= 0.5)
  
  # (3) ADD LABELS
  cc.cxc.ack <- read.csv("input/cc_cxc_ack.csv")
  data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
  rm(cc.cxc.ack)
  
  # (4) COUNT HOW MANY CC AND CXC MOTIFS you are querying from FOR NORMALIZATION
  how.many <- data %>% filter(pct_ortho >= 0.5)
  how.many <- how.many %>% select(class) %>% count(class) %>% 
    filter(class == "cc" | class == "cxc")
  print(how.many) # 288 CC, 169 CXC
  
  # (5) COUNT HOW MANY MOTIFS CONTAIN AROMATIC RESIDUES
  tyr <- subset(data, grepl("Y", data$motif))
  tyr$aa <- c("tyr")
  phe <- subset(data, grepl("F", data$motif))
  phe$aa <- c("phe")
  pro <- subset(data, grepl("P", data$motif))
  pro$aa <- c("pro")
  trp <- subset(data, grepl("W", data$motif))
  trp$aa <- c("trp")
  
  # (7) BIND CHARGED MOTIFS TOGETHER & ADD PLUS MINUS
  arom <- bind_rows(tyr, phe, pro, trp)
  #rm(asp, glu, arg, lys, his)
  
  
  # remove residues counted twice
  arom <- arom %>% select(motif, protein, class)
  arom <- arom %>% distinct()
  
  # (8) NORMALIZE BY NUMBER OF RESIDUES QUERIED WITHIN EACH FAMILY
  class <- arom %>% count(class)
  arom <- arom %>%
    count(class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/288,
      class == "cxc" ~ n/169
    ))
  
  # add dummy column
  arom$dummy <- c("dummy")
  
  # (9) BAR PLOT - CHARGE
  arom %>%
    filter(class != "cx3c" & class != "xc" & class != "ack") %>%
    ggplot(aes(dummy, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.6) +
    theme_minimal() +
    facet_grid(. ~ class)
  