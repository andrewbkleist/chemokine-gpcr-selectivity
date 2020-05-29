# Name:     7_tier_and_glay_pie.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Supp. Figure 12A,B

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
  
  # override 3x32 as "honorary" Tier2
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
  
  # remove "_" notation
  rin.sum <- rin.sum %>% mutate(target_gnccn = case_when(
    target_gnccn == "1x22_" ~ "1x22",
    target_gnccn == "1x23_" ~ "1x23",
    target_gnccn == "1x24_" ~ "1x24",
    target_gnccn == "1x28_" ~ "1x28",
    target_gnccn != "1x22_" & target_gnccn != "1x23_" & target_gnccn != "1x24_" & target_gnccn != "1x28_" ~ target_gnccn 
  ))
  
  # (4) MAP GLAY TO CONTACTS ---------------------------------------------------
  
  # import glay clusters
  glay <- read_csv("output/old/CONTACT_TIER_GLAY_5UIW.csv") 
    # has old Tier designtaions - must remove...
  glay <- glay %>% select(-c(ck_tier:ckr_tier_adj))

  # add 1x22, 1x23
  glay <- glay %>% mutate(target_gnccn = case_when(
    target_gnccn == "NTr.Cys" ~ "1x22",
    target_gnccn == "NTr.Cp1" ~ "1x23",
    target_gnccn != "1x22" & target_gnccn != "1x23" ~ target_gnccn 
  ))
  
  # combine glay clusters with RIN/tier info
  rin.sum$glay_ck <- glay$glayck[match(unlist(rin.sum$source_gnccn), glay$source_gnccn)]
  rin.sum$glay_ckr <- glay$glayckr[match(unlist(rin.sum$target_gnccn), glay$target_gnccn)]
  
  # remove
  rm(glay, lig, rec)
  
  # (5) MAP GLAY TO CONTACTS ---------------------------------------------------
  rin.sum <- rin.sum %>% select(tier_ck, tier_ckr, glay_ck, glay_ckr)
  
  # get numbers of contacts by tier per each community cluster
  ck <- rin.sum %>% count(glay_ck, tier_ck)
  ck$type <- c("ck")
  colnames(ck) <- c("module", "tier", "n", "type")
  ckr <- rin.sum %>% count(glay_ckr, tier_ckr)
  ckr$type <- c("ckr")
  colnames(ckr) <- c("module", "tier", "n", "type")
  master <- bind_rows(ck, ckr)
  
  # remove
  rm(ck, ckr, rin.sum)
  
  # change designations
  master <- master %>% mutate(tier = case_when(
    tier == "NTc" ~ "unstr",
    tier == "NTr" ~ "unstr",
    tier == "ECL3" ~ "other",
    tier != "NTc" & tier != "NTr" & tier != "ECL3" ~ tier 
  ))
  
  # rearrange data for pie chart and get normalized percentages by number of contacts
  master2 <- master %>% group_by(module, type) %>% 
    mutate(mod_sum =(
      sum(n)
    )) %>%
    ungroup()
  master2 <- master2 %>% mutate(mod_norm = n / mod_sum)
  
  # order community clusters
  order <- c(1:8)
  master2$module <- factor(master2$module, levels = order)
  
  master2 %>%
    ggplot(aes(x = "", y= mod_norm, fill = tier)) +
    geom_bar(width = 1, size = 0.1, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_manual(values=c(ECL2 = "purple4", 
                               other = "white", 
                               tier1 = "gray20",
                               tier2 = "dodgerblue4", 
                               tier3 = "mediumslateblue",
                               unstr = "purple4",
                               low_ortho = "grey85")) +
    facet_grid(module ~ type)
  