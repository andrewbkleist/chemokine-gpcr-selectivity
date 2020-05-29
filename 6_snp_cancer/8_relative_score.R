# Name:     8_relative_score.R
# Updated:  20191102
# Author:   Greg Slodkowicz / Andrew Kleist
# Figure:   Figure 6I,J

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################
  
  library(tidyverse)
  library(ggrepel)
  library(ggforce)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")

##### 1: CHEMOKINE ANALYSIS ####################################################
  
  # (1) IMPORT DATA ------------------------------------------------------------
  data <- read_csv("output/CK_CANCER_SNP_TIER_OUTPUT.csv")
  data <- data %>% mutate(selectivity = case_when(
    tier == "motif" ~ "yes",
    tier == "fragment" ~ "no",
    tier == "tier1" ~ "yes",
    tier == "tier2" ~ "yes",
    tier == "tier3" ~ "yes",
    tier == "low_ortho" ~ "no")) %>%
    filter(count < 30000) 
  # an allele with a count of 30,000 is 10% of 300,000 alleles in GNOMAD
  
  # (2) DEFINE FUNCTIONS
  GetPositionsA <- function(DF, TIER, TYPE){
    data <- DF %>% 
      filter(type == TYPE) %>%
      filter(tier == TIER) %>%
      filter(pctgaps < 0.6) %>%
      dplyr::summarise(n=dplyr::n(), count = sum(count))
    data$tier <- c(TIER)
    data$type <- c(TYPE)
    data <- data %>% mutate(rel = count / n)
    return(data)
  }
 
  GetPositionsB <- function(DF, TIER, TYPE){
    data <- DF %>% 
      filter(type == TYPE) %>%
      filter(tier == TIER) %>%
      #filter(pctgaps < 0.6) %>%
      dplyr::summarise(n=dplyr::n(), count = sum(count))
    data$tier <- c(TIER)
    data$type <- c(TYPE)
    data <- data %>% mutate(rel = count / n)
    return(data)
  }
  
  # (3) GET SNPs
  
  # get SNP structured
  tier1 <- GetPositionsA(data, "tier1", "snp_count")
  tier2 <- GetPositionsA(data, "tier2", "snp_count")
  tier3 <- GetPositionsA(data, "tier3", "snp_count")
  nontier <- GetPositionsA(data, "low_ortho", "snp_count")
  ck.snp.str <- bind_rows(tier1, tier2, tier3) %>%
    mutate(score = (rel/nontier$rel), logscore = log2(rel/nontier$rel))
  rm(tier1, tier2, tier3, nontier)
  
  # get SNP unstructured
  motif <- GetPositionsB(data, "motif", "snp_count")
  fragment <- GetPositionsB(data, "fragment", "snp_count")
  ck.snp.unstr <- motif %>%
    mutate(score = (rel/fragment$rel), logscore = log2(rel/fragment$rel))
  rm(motif, fragment)
  
  # combine
  snp <- bind_rows(ck.snp.str, ck.snp.unstr)
  rm(ck.snp.str, ck.snp.unstr)
  
  # (4) GET CANCER
  tier1 <- GetPositionsA(data, "tier1", "cancer_count")
  tier2 <- GetPositionsA(data, "tier2", "cancer_count")
  tier3 <- GetPositionsA(data, "tier3", "cancer_count")
  nontier <- GetPositionsA(data, "low_ortho", "cancer_count")
  ck.can.str <- bind_rows(tier1, tier2, tier3) %>%
    mutate(score = (rel/nontier$rel), logscore = log2(rel/nontier$rel))
  rm(tier1, tier2, tier3, nontier)
  
  # get SNP unstructured
  motif <- GetPositionsB(data, "motif", "cancer_count")
  fragment <- GetPositionsB(data, "fragment", "cancer_count")
  ck.can.unstr <- motif %>%
    mutate(score = (rel/fragment$rel), logscore = log2(rel/fragment$rel))
  rm(motif, fragment)
  
  # combine
  can <- bind_rows(ck.can.str, ck.can.unstr)
  rm(ck.can.str, ck.can.unstr)

  
  # (5) ADD COMBINED NUMBERS
  ck.snp <- log2( (27522/1700) / (62222/1189) )
  ck.can <- log2( (311/1700) / (269/1194) )
  df <- data.frame(n = NA, count = c(NA, NA) , tier = c("all", "all"), 
         type = c("snp_count", "cancer_count"), rel = c(NA, NA), 
         score = c(NA, NA), logscore = c(ck.snp, ck.can) )
    
  # (6) COMBINE ALL & GRAPH
  all.ck <- bind_rows(snp, can, df)
  rm(df, ck.snp, ck.can)
  
  order <- c("tier1", "tier2", "tier3", "motif", "all")
  all.ck$tier <- factor(all.ck$tier, levels = order)
  order2 <- c("snp_count", "cancer_count")
  all.ck$type <- factor(all.ck$type, levels = (order2))
  rm(snp, can)
  
  all.ck %>%
    ggplot(aes(tier, logscore, group=type, fill=type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylim(-5, 0.6) +
    theme_minimal()
  
  
  
##### 2: RECEPTOR ANALYSIS #####################################################
  
  # (1) IMPORT DATA ------------------------------------------------------------
  data <- read.csv("output/CKR_CANCER_SNP_TIER_OUTPUT.csv")
  data <- data %>% mutate(selectivity = case_when(
    tier == "motif" ~ "yes",
    tier == "fragment" ~ "no",
    tier == "tier1" ~ "yes",
    tier == "tier2" ~ "yes",
    tier == "tier3" ~ "yes",
    tier == "non_tier" ~ "no")) %>%
    filter(count < 30000)
  # an allele with a count of 30,000 is 10% of 300,000 alleles in GNOMAD
  
  # (2) FUNCTIONS (above)
  
  # (3) GET SNPs
  # get SNP structured
  tier1 <- GetPositionsA(data, "tier1", "snp_count")
  tier2 <- GetPositionsA(data, "tier2", "snp_count")
  tier3 <- GetPositionsA(data, "tier3", "snp_count")
  nontier <- GetPositionsA(data, "non_tier", "snp_count")
  ck.snp.str <- bind_rows(tier1, tier2, tier3) %>%
    mutate(score = (rel/nontier$rel), logscore = log2(rel/nontier$rel))
  rm(tier1, tier2, tier3, nontier)
  
  # get SNP unstructured
  motif <- GetPositionsB(data, "motif", "snp_count")
  fragment <- GetPositionsB(data, "fragment", "snp_count")
  ck.snp.unstr <- motif %>%
    mutate(score = (rel/fragment$rel), logscore = log2(rel/fragment$rel))
  rm(motif, fragment)
  
  # combine
  snp <- bind_rows(ck.snp.str, ck.snp.unstr)
  rm(ck.snp.str, ck.snp.unstr)
  
  # (4) GET CANCER
  tier1 <- GetPositionsA(data, "tier1", "cancer_count")
  tier2 <- GetPositionsA(data, "tier2", "cancer_count")
  tier3 <- GetPositionsA(data, "tier3", "cancer_count")
  nontier <- GetPositionsA(data, "non_tier", "cancer_count")
  ck.can.str <- bind_rows(tier1, tier2, tier3) %>%
    mutate(score = (rel/nontier$rel), logscore = log2(rel/nontier$rel))
  rm(tier1, tier2, tier3, nontier)
  
  # get SNP unstructured
  motif <- GetPositionsB(data, "motif", "cancer_count")
  fragment <- GetPositionsB(data, "fragment", "cancer_count")
  ck.can.unstr <- motif %>%
    mutate(score = (rel/fragment$rel), logscore = log2(rel/fragment$rel))
  rm(motif, fragment)
  
  # combine
  can <- bind_rows(ck.can.str, ck.can.unstr)
  rm(ck.can.str, ck.can.unstr)
  
  # (5) ADD COMBINED NUMBERS
  ckr.snp <- log2( (10014/1885) / (78867/3051) )
  ckr.can <- log2( (378/1886) / (499/3052) )
  df <- data.frame(n = NA, count = c(NA, NA) , tier = c("all", "all"), 
                   type = c("snp_count", "cancer_count"), rel = c(NA, NA), 
                   score = c(NA, NA), logscore = c(ckr.snp, ckr.can) )
  
  # (6) COMBINE ALL & GRAPH
  all.ckr <- bind_rows(snp, can, df)
  rm(df, ckr.snp, ckr.can)
  
  order <- c("tier1", "tier2", "tier3", "motif", "all")
  all.ckr$tier <- factor(all.ckr$tier, levels = order)
  order2 <- c("snp_count", "cancer_count")
  all.ckr$type <- factor(all.ckr$type, levels = (order2))
  rm(snp, can)
  
  all.ckr %>%
    ggplot(aes(tier, logscore, group=type, fill=type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylim(-5, 0.6) +
    theme_minimal()
  
  