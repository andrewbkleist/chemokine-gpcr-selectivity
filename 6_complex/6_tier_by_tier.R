# Name:     6_tier_by_tier.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Figure 

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ggrepel)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_complex/")
  
##### 1: MAP TIERS TO ALL-BY-ALL INTERACTION MATRIX ############################
  
  # (1) CONTACTS ---------------------------------------------------------------
  # import contacts, select interface only, select Xray only
  rin <- read_csv("input/RIN.csv") %>%
    mutate(type = case_when(
      Chain1 == "A" & Chain2 == "A" ~ "ck",
      Chain1 == "B" & Chain2 == "B" ~ "ckr",
      Chain1 == "A" & Chain2 == "B" ~ "inter"
    )) %>% filter(type =="inter") %>%
    filter(class == "xray") %>%
    filter(file == "5uiw")
  
  # simplify DF
  rin <- rin %>% select(-PDB, -Chain1, -Chain2, -SS1, -SS2, 
                        -Number.of.atomic.contacts, -class, -type)
  
  # count occurances of individual contacts
  rin.sum <- rin %>% select(source_gnccn, target_gnccn) %>% count(source_gnccn, target_gnccn)
  rm(rin)
  
  # (2) TIERS ------------------------------------------------------------------
  lig <- read.csv("../2_ck_core_tier/output/MASTER_CONSERVATION_WITH_TIERS.csv") 
  lig <- lig %>% filter(ck == "ccl5") %>% filter(!is.na(tier))
  
  rec <- read.csv("../4_ckr_core_tier/output/MASTER_CONSERVATION_WITH_TIERS.csv")
  rec <- rec %>% filter(ckr == "ccr5") %>% filter(!is.na(tier)) %>% 
    filter(dom != "NT" & dom != "CT" & dom != "ICL2"  & dom != "ICL3") %>%
    filter(ortho_cons != 0)
  
  # (3) MAP TIERS TO CONTACTS --------------------------------------------------
  # remove "gn" from GPCRdb names
  rec <- rec %>% separate(gn, into = c("gn1", "gn"), sep = "gn") %>% select(-gn1)
  
  # equate "1x24" and "1x24_" nomenclature
  rin.sum <- rin.sum %>% mutate(target_gnccn = case_when(
    target_gnccn == "NTr.Cys" ~ "1x22_",
    target_gnccn == "NTr.Cp1" ~ "1x23_",
    target_gnccn == "1x24" ~ "1x24_",
    target_gnccn == "1x25" ~ "1x25_",
    target_gnccn == "1x26" ~ "1x26_",
    target_gnccn == "1x27" ~ "1x27_",
    target_gnccn == "1x28" ~ "1x28_",
    target_gnccn == "1x29" ~ "1x29_",
    target_gnccn != "1x24" & target_gnccn != "NTr.Cp1" & target_gnccn != "NTr.Cys" ~ target_gnccn
  ))
  
  # match res and fill Tier for chemokine and receptor
  rin.sum$tier_ck <- lig$tier[match(unlist(rin.sum$source_gnccn), lig$ccn)]
  rin.sum$tier_ckr <- rec$tier[match(unlist(rin.sum$target_gnccn), rec$gn)]
  
  # override 3x32 as Tier2 (see methods)
  rin.sum$tier_ckr <- as.character(rin.sum$tier_ckr)
  rin.sum <- rin.sum %>% mutate(tier_ckr = case_when(
    target_gnccn == "3x32" ~ "tier2",
    target_gnccn != "3x32" ~ tier_ckr
  ))
  
  # add designations for unstructured regions
  rin.sum <- rin.sum %>% mutate(tier_ckr = case_when(
    !is.na((tier_ckr)) ~ as.character(tier_ckr),
    grepl("NTr", rin.sum$target_gnccn) ~ "NTr",
    grepl("ECL2", rin.sum$target_gnccn) ~ "ECL2",
    grepl("ECL3", rin.sum$target_gnccn) ~ "ECL3"
  )) %>%
    mutate(tier_ck = case_when(
      !is.na((tier_ck)) ~ as.character(tier_ck),
      grepl("NTc", rin.sum$source_gnccn) ~ "NTc"
    ))

  # (4) SUMMARY PLOT -----------------------------------------------------------
  # count tier-to-tier contacts
  tier.sum <- rin.sum %>% count(tier_ck, tier_ckr)
  
  # order
  order.ck <- c("tier1", "tier2", "tier3", "low_ortho", "NTc")
  order.ckr <- c("tier1", "tier2", "tier3", "low_ortho", "ECL3", "NTr", "ECL2")
  tier.sum$tier_ck <- factor(tier.sum$tier_ck, levels = rev(order.ck))
  tier.sum$tier_ckr <- factor(tier.sum$tier_ckr, levels = order.ckr)
  
  # TILE PLOT
  tier.sum %>%
    ggplot(aes(tier_ckr, tier_ck, fill = n)) +
    geom_tile() +
    scale_fill_gradient(low="grey90", high="black") +
    geom_text(aes(tier_ckr, tier_ck, label = n), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  
  
  