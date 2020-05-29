# Name:     5_arg_vs_asp.R
# Updated:  20191031
# User:     Andrew Kleist
# Figure:   Figure 3D; Supp. Figure 6A

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/3_ck_motif/")

##### (1) GET CHARGE DISTRIBUTION AND PROLINE FOR ALL MERS #########################

  # (1) PREP DATA FRAME
  # import
  data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")
  
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
  print(how.many) # 923 CC, 638 CXC
  
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
  
  # (6) COUNT HOW MANY MOTIFS CONTAIN PROLINES
  pro <- subset(data, grepl("P", data$motif))
  pro$aa <- c("pro")
  
  # (7) BIND CHARGED MOTIFS TOGETHER & ADD PLUS MINUS
  charged <- bind_rows(asp, glu, lys, arg, his)
  rm(asp, glu, arg, lys, his)
  
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
  charged <- charged %>% filter(pct_ortho >= 0.5)
  class <- charged %>% count(class)
  charged <- charged %>%
    count(charge, class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/923,
      class == "cxc" ~ n/638
    ))
  
  # (9) BAR PLOT - CHARGE
  charged %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(charge, nl)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(0, 0.4) +
    facet_grid(. ~ class)
  
  # (10) PROLINE NORMALIZATION
  pro <- pro %>% filter(pct_ortho >= 0.5)
  class <- pro %>% count(class)
  pro <- pro %>%
    count(class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/923,
      class == "cxc" ~ n/638
    ))
  
  # (11) BAR PLOT - PROLINE
  pro %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(class, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.4) +
    theme_minimal() 
  
##### (2) GET CHARGE DISTRIBUTION AND PROLINE FOR ALL 2 MERS INIDIVIDUALLY #####
  
  # (1) PREP DATA FRAME
  # import
  data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")
  
  # select only CONSERVED MOTIFS
  data <- data %>% filter(pct_ortho >= 0.5)
  
  # SELECT ONLY 2-MER AND NOT MASKED MOTIFS
  data <- data %>% filter(mer ==  "mer2")
  data <- data %>% filter(mask == "none")
  
  # (3) ADD LABELS
  cc.cxc.ack <- read.csv("input/cc_cxc_ack.csv")
  data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
  rm(cc.cxc.ack)
  
  # (4) COUNT HOW MANY CC AND CXC MOTIFS you are querying from FOR NORMALIZATION
  how.many <- data %>% filter(pct_ortho >= 0.5)
  how.many <- how.many %>% select(class) %>% count(class) %>% 
    filter(class == "cc" | class == "cxc")
  print(how.many) # 181 CC, 122 CXC
  
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
  
  # (6) COUNT HOW MANY MOTIFS CONTAIN PROLINES
  pro <- subset(data, grepl("P", data$motif))
  pro$aa <- c("pro")
  
  # (7) BIND CHARGED MOTIFS TOGETHER & ADD PLUS MINUS
  charged <- bind_rows(asp, glu, lys, arg, his)
  rm(asp, glu, arg, lys, his)
  
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
  charged <- charged %>% filter(pct_ortho >= 0.5)
  class <- charged %>% count(class)
  charged <- charged %>%
    count(charge, class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/181,
      class == "cxc" ~ n/122
    ))
  
  # (9) BAR PLOT - CHARGE
  charged %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(charge, nl)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(0, 0.6) +
    facet_grid(. ~ class)
  
  # (10) PROLINE NORMALIZATION
  pro <- pro %>% filter(pct_ortho >= 0.5)
  class <- pro %>% count(class)
  pro <- pro %>%
    count(class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/181,
      class == "cxc" ~ n/122
    ))
  
  # (11) BAR PLOT - PROLINE
  pro %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(class, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.6) +
    theme_minimal() 
  
##### (3) GET CHARGE DISTRIBUTION AND PROLINE FOR ALL 3 MERS INIDIVIDUALLY #####
  
  # (1) PREP DATA FRAME
  # import
  data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")
  
  # select only CONSERVED MOTIFS
  data <- data %>% filter(pct_ortho >= 0.5)
  
  # SELECT ONLY 2-MER AND NOT MASKED MOTIFS
  data <- data %>% filter(mer ==  "mer3")
  data <- data %>% filter(mask == "none")
  
  # (3) ADD LABELS
  cc.cxc.ack <- read.csv("input/cc_cxc_ack.csv")
  data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
  rm(cc.cxc.ack)
  
  # (4) COUNT HOW MANY CC AND CXC MOTIFS you are querying from FOR NORMALIZATION
  how.many <- data %>% filter(pct_ortho >= 0.5)
  how.many <- how.many %>% select(class) %>% count(class) %>% 
    filter(class == "cc" | class == "cxc")
  print(how.many) # 123 CC, 94 CXC
  
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
  
  # (6) COUNT HOW MANY MOTIFS CONTAIN PROLINES
  pro <- subset(data, grepl("P", data$motif))
  pro$aa <- c("pro")
  
  # (7) BIND CHARGED MOTIFS TOGETHER & ADD PLUS MINUS
  charged <- bind_rows(asp, glu, lys, arg, his)
  rm(asp, glu, arg, lys, his)
  
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
  charged <- charged %>% filter(pct_ortho >= 0.5)
  class <- charged %>% count(class)
  charged <- charged %>%
    count(charge, class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/123,
      class == "cxc" ~ n/94
    ))
  
  # (9) BAR PLOT - CHARGE
  charged %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(charge, nl)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(0, 0.6) +
    facet_grid(. ~ class)
  
  # (10) PROLINE NORMALIZATION
  pro <- pro %>% filter(pct_ortho >= 0.5)
  class <- pro %>% count(class)
  pro <- pro %>%
    count(class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/123,
      class == "cxc" ~ n/94
    ))
  
  # (11) BAR PLOT - PROLINE
  pro %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(class, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.6) +
    theme_minimal() 

    
##### (4) GET CHARGE DISTRIBUTION AND PROLINE FOR ALL 4 MERS INIDIVIDUALLY #####
  
  # (1) PREP DATA FRAME
  # import
  data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")
  
  # select only CONSERVED MOTIFS
  data <- data %>% filter(pct_ortho >= 0.5)
  
  # SELECT ONLY 2-MER AND NOT MASKED MOTIFS
  data <- data %>% filter(mer ==  "mer4")
  data <- data %>% filter(mask == "none")
  
  # (3) ADD LABELS
  cc.cxc.ack <- read.csv("input/cc_cxc_ack.csv")
  data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
  rm(cc.cxc.ack)
  
  # (4) COUNT HOW MANY CC AND CXC MOTIFS you are querying from FOR NORMALIZATION
  how.many <- data %>% filter(pct_ortho >= 0.5)
  how.many <- how.many %>% select(class) %>% count(class) %>% 
    filter(class == "cc" | class == "cxc")
  print(how.many) # 90 CC, 66 CXC
  
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
  
  # (6) COUNT HOW MANY MOTIFS CONTAIN PROLINES
  pro <- subset(data, grepl("P", data$motif))
  pro$aa <- c("pro")
  
  # (7) BIND CHARGED MOTIFS TOGETHER & ADD PLUS MINUS
  charged <- bind_rows(asp, glu, lys, arg, his)
  rm(asp, glu, arg, lys, his)
  
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
  charged <- charged %>% filter(pct_ortho >= 0.5)
  class <- charged %>% count(class)
  charged <- charged %>%
    count(charge, class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/90,
      class == "cxc" ~ n/66
    ))
  
  # (9) BAR PLOT - CHARGE
  charged %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(charge, nl)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(0, 0.6) +
    facet_grid(. ~ class)
  
  # (10) PROLINE NORMALIZATION
  pro <- pro %>% filter(pct_ortho >= 0.5)
  class <- pro %>% count(class)
  pro <- pro %>%
    count(class) %>%
    mutate(nl = case_when(
      class == "cc" ~ n/90,
      class == "cxc" ~ n/66
    ))
  
  # (11) BAR PLOT - PROLINE
  pro %>%
    filter(class != "cx3c" & class != "xc") %>%
    ggplot(aes(class, nl)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.6) +
    theme_minimal() 