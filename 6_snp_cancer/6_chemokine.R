# Name:     6_chemokine.R
# Updated:  20191219
# Author:   Greg Slodkowicz / Andrew Kleist
# Figure:   Figure 6I,J, Supp. Figure 13

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################
  
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")

##### 1: CHEMOKINE ANALYSIS - SNP ##############################################

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
    filter(tier == "tier1" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(410, 40322), p=c(252, 880)/(252 + 880))
  tier1 <- ChiToGraph(tier1, "tier1", "low_ortho")
    # test <- tier1 %>% spread(EO, value)
    # test <- test %>% mutate(diff = O - E, chi2 = ((O-E)^2)/E )
  
  order = c("tier1", "low_ortho")
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
    filter(tier == "tier2" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(9252, 40322), p=c(829, 880)/(829 + 880))
  tier2 <- ChiToGraph(tier2, "tier2", "low_ortho")
  
  order = c("tier2", "low_ortho")
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
    filter(tier == "tier3" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(14196, 40322), p=c(431, 880)/(431 + 880))
  tier3 <- ChiToGraph(tier3, "tier3", "low_ortho")
  
  order = c("tier3", "low_ortho")
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
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(23858, 40322), p=c(1512, 880)/(1512 + 880))
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
  
  # (4) NTERM
  nt <- data %>% filter(type == "snp_count") %>%
    filter(tier == "motif" | tier == "fragment") %>%
    #filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(3664, 21900), p=c(188, 309)/(188 + 309))
  nt <- ChiToGraph(nt, "motif", "fragment")
  
  order = c("motif", "fragment")
  order2 = c("E", "O")
  nt$tier <- factor(nt$tier, levels = order)
  nt$EO <- factor(nt$EO, levels = order2)
  nt %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("purple", "grey70")) +
    theme_minimal()
  
  # (5) GROUPED
  alltier <- data %>% filter(type == "snp_count") %>%
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) 
  motif.frag <- data %>% filter(type == "snp_count") %>%
    filter(tier == "motif" | tier == "fragment")
  all <- bind_rows(alltier, motif.frag)
  all <- all %>% group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(27522, 62222), p=c(1700, 1189)/(1700 + 1189))
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
  
##### 2: CHEMOKINE ANALYSIS - CANCER  ##########################################
  
  # (1) IMPORT DATA ------------------------------------------------------------
  data <- read_csv("output/CK_CANCER_SNP_TIER_OUTPUT.csv")
  data <- data %>% mutate(selectivity = case_when(
    tier == "motif" ~ "yes",
    tier == "fragment" ~ "no",
    tier == "tier1" ~ "yes",
    tier == "tier2" ~ "yes",
    tier == "tier3" ~ "yes",
    tier == "low_ortho" ~ "no"))

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
    filter(tier == "tier1" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(36, 205), p=c(252, 883)/(252 + 883))
  tier1 <- ChiToGraph(tier1, "tier1", "low_ortho")
  # test <- tier1 %>% spread(EO, value)
  # test <- test %>% mutate(diff = O - E, chi2 = ((O-E)^2)/E )
  
  order = c("tier1", "low_ortho")
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
    filter(tier == "tier2" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(152, 205), p=c(829, 883)/(829 + 883))
  tier2 <- ChiToGraph(tier2, "tier2", "low_ortho")
  
  order = c("tier2", "low_ortho")
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
    filter(tier == "tier3" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(96, 205), p=c(431, 883)/(431 + 883))
  tier3 <- ChiToGraph(tier3, "tier3", "low_ortho")
  
  order = c("tier3", "low_ortho")
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
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count))
  chisq.test(c(23858, 40322), p=c(1512, 880)/(1512 + 880))
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
  
  # (4) NTERM
  nt <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "motif" | tier == "fragment") %>%
    #filter(pctgaps < 0.6) %>%
    #filter(ortho_cons > 0.8) %>%
    group_by(tier) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(27, 64), p=c(188, 311)/(188 + 311))
  nt <- ChiToGraph(nt, "motif", "fragment")
  
  order = c("motif", "fragment")
  order2 = c("E", "O")
  nt$tier <- factor(nt$tier, levels = order)
  nt$EO <- factor(nt$EO, levels = order2)
  nt %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("purple", "grey70")) +
    theme_minimal()
  
  # (5) GROUPED
  alltier <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "tier1" | tier == "tier2" | tier == "tier3" | tier == "low_ortho") %>%
    filter(pctgaps < 0.6) 
  motif.frag <- data %>% filter(type == "cancer_count") %>%
    filter(tier == "motif" | tier == "fragment")
  all <- bind_rows(alltier, motif.frag)
  all <- all %>% group_by(selectivity) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(count)) %>% ungroup()
  chisq.test(c(311, 269), p=c(1700, 1194)/(1700 + 1194))
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
  