# Name:     7_receptor.R
# Updated:  20191219
# Author:   Greg Slodkowicz / Andrew Kleist
# Figure:   Figure 6I,J, Supp. Figure 13

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################
  
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")

##### 1: RECEPTOR ANALYSIS - SNP ###############################################
# Comparisons of SNPs/cancer muts from Tier versus non...

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
  
  # (2) DEFINE FUNCTION FOR EXPECTED VS. ACTUAL GRAPH
  ChiToGraph <- function(DF, TIER, NON){
    data <- DF %>% gather(key, value, 2:3)
    data <- data %>% dplyr::group_by(key) %>% dplyr::mutate(total = sum(value)) %>% dplyr::ungroup()
    data <- data %>% mutate(fraction = value / total)
    data <- data.frame("EO" = c("O", "E", "O", "E"), 
                       "tier" = c(TIER, TIER, NON, NON),
                       "value" = c(data$value[4], data$fraction[2]*data$total[3], data$value[3], data$fraction[1]*data$total[3]) )
    return(data)
    rm(data)
  }
  
  # (3) CHI SQUARE TEST
  # TIER 1
  tier1 <- data %>% filter(type == "snp_count") %>%
    filter(tier == "tier1" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(669, 57323), p=c(303, 2442)/(303 + 2442))
  tier1 <- ChiToGraph(tier1, "tier1", "non_tier")
  
  order = c("tier1", "non_tier")
  order2 = c("E", "O")
  tier1$tier <- factor(tier1$tier, levels = order)
  tier1$EO <- factor(tier1$EO, levels = order2)
  tier1 %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # TIER 2
  tier2 <- data %>% filter(type == "snp_count") %>%
    filter(tier == "tier2" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(5367, 57323), p=c(575, 2442)/(575 + 2442))
  tier2 <- ChiToGraph(tier2, "tier2", "non_tier")
  
  order = c("tier2", "non_tier")
  order2 = c("E", "O")
  tier2$tier <- factor(tier2$tier, levels = order)
  tier2$EO <- factor(tier2$EO, levels = order2)
  tier2 %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("steelblue3", "grey70")) +
    theme_minimal()
  
  # TIER 3
  tier3 <- data %>% filter(type == "snp_count") %>%
    filter(tier == "tier3" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(3749, 57323), p=c(879, 2442)/(879 + 2442))
  tier3 <- ChiToGraph(tier3, "tier3", "non_tier")
  
  order = c("tier3", "non_tier")
  order2 = c("E", "O")
  tier3$tier <- factor(tier3$tier, levels = order)
  tier3$EO <- factor(tier3$EO, levels = order2)
  tier3 %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("thistle3", "grey70")) +
    theme_minimal()

  # ALL TIERS
  alltier <- data %>% filter(type == "snp_count") %>%
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(9785, 57323), p=c(1756, 2442)/(1756 + 2442))
  alltier <- ChiToGraph(alltier, "yes", "no")
  
  order = c("yes", "no")
  order2 = c("E", "O")
  alltier$tier <- factor(alltier$tier, levels = order)
  alltier$EO <- factor(alltier$EO, levels = order2)
  alltier %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # (4) NTERM & ECL2
  unst <- data %>% filter(type == "snp_count") %>%
    filter(tier == "motif" | tier == "fragment") %>%
    #filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(231, 21547), p=c(129, 614)/(129 + 614))
  unst <- ChiToGraph(unst, "yes", "no")
  
  order = c("yes", "no")
  order2 = c("E", "O")
  unst$tier <- factor(unst$tier, levels = order)
  unst$EO <- factor(unst$EO, levels = order2)
  unst %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("purple", "grey70")) +
    theme_minimal()
  
  # (5) GROUPED
  alltier <- data %>% filter(type == "snp_count") %>%
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) 
  motif.frag <- data %>% filter(type == "snp_count") %>%
    filter(tier == "motif" | tier == "fragment")
  all <- bind_rows(alltier, motif.frag)
  all <- all %>% group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(10016, 78870), p=c(1885, 3056)/(1885 + 3056))
  all <- ChiToGraph(all, "yes", "no")
  
  order = c("yes", "no")
  order2 = c("E", "O")
  all$tier <- factor(all$tier, levels = order)
  all$EO <- factor(all$EO, levels = order2)
  all %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  
##### 2: RECEPTOR ANALYSIS - CANCER ############################################

  # (1) IMPORT DATA ------------------------------------------------------------
  data <- read.csv("output/CKR_CANCER_SNP_TIER_OUTPUT.csv")
  data <- data %>% mutate(selectivity = case_when(
    tier == "motif" ~ "yes",
    tier == "fragment" ~ "no",
    tier == "tier1" ~ "yes",
    tier == "tier2" ~ "yes",
    tier == "tier3" ~ "yes",
    tier == "non_tier" ~ "no")) 

  # (2) DEFINE FUNCTION FOR EXPECTED VS. ACTUAL GRAPH
  ChiToGraph <- function(DF, TIER, NON){
    data <- DF %>% gather(key, value, 2:3)
    data <- data %>% dplyr::group_by(key) %>% dplyr::mutate(total = sum(value)) %>% dplyr::ungroup()
    data <- data %>% mutate(fraction = value / total)
    data <- data.frame("EO" = c("O", "E", "O", "E"), 
                       "tier" = c(TIER, TIER, NON, NON),
                       "value" = c(data$value[4], data$fraction[2]*data$total[3], data$value[3], data$fraction[1]*data$total[3]) )
    return(data)
    rm(data)
  }
  
  # (3) CHI SQUARE TEST
  # TIER 1
  tier1 <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "tier1" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(49, 364), p=c(303, 2443)/(303 + 2443))
  tier1 <- ChiToGraph(tier1, "tier1", "non_tier")
  
  order = c("tier1", "non_tier")
  order2 = c("E", "O")
  tier1$tier <- factor(tier1$tier, levels = order)
  tier1$EO <- factor(tier1$EO, levels = order2)
  tier1 %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # TIER 2
  tier2 <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "tier2" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(123, 364), p=c(575, 2443)/(575 + 2443))
  tier2 <- ChiToGraph(tier2, "tier2", "non_tier")
  
  order = c("tier2", "non_tier")
  order2 = c("E", "O")
  tier2$tier <- factor(tier2$tier, levels = order)
  tier2$EO <- factor(tier2$EO, levels = order2)
  tier2 %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("steelblue3", "grey70")) +
    theme_minimal()
  
  # TIER 3
  tier3 <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "tier3" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(184, 364), p=c(879, 2443)/(879 + 2443))
  tier3 <- ChiToGraph(tier3, "tier3", "non_tier")
  
  order = c("tier3", "non_tier")
  order2 = c("E", "O")
  tier3$tier <- factor(tier3$tier, levels = order)
  tier3$EO <- factor(tier3$EO, levels = order2)
  tier3 %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("thistle3", "grey70")) +
    theme_minimal()
  
  # ALL TIERS
  alltier <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(9785, 57323), p=c(1756, 2442)/(1756 + 2442))
  alltier <- ChiToGraph(alltier, "yes", "no")
  
  order = c("yes", "no")
  order2 = c("E", "O")
  alltier$tier <- factor(alltier$tier, levels = order)
  alltier$EO <- factor(alltier$EO, levels = order2)
  alltier %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # (4) NTERM & ECL2
  unst <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "motif" | tier == "fragment") %>%
    #filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(22, 136), p=c(129, 614)/(129 + 614))
  unst <- ChiToGraph(unst, "yes", "no")
  
  order = c("yes", "no")
  order2 = c("E", "O")
  unst$tier <- factor(unst$tier, levels = order)
  unst$EO <- factor(unst$EO, levels = order2)
  unst %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("purple", "grey70")) +
    theme_minimal()
  
  # (5) GROUPED
  alltier <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "non_tier") %>%
    filter(pctgaps < 0.6) 
  motif.frag <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "motif" | tier == "fragment")
  all <- bind_rows(alltier, motif.frag)
  all <- all %>% group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(378, 500), p=c(1886, 3057)/(1886 + 3057))
  all <- ChiToGraph(all, "yes", "no")
  
  order = c("yes", "no")
  order2 = c("E", "O")
  all$tier <- factor(all$tier, levels = order)
  all$EO <- factor(all$EO, levels = order2)
  all %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
