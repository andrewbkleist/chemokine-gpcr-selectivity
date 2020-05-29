# Name:     2_unstr_str_cyto.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Figure 6E

# This script makes the input for Figure 6E, which was loaded 
# and plotted in cytoscape

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_complex/")

##### STRUCTURED AND UNSTRUCTURED CYTOSCAPE ####################################
  
  # (1) IMPORT CONTACTS
  
  # import contacts, select interface only, select Xray only
  rin <- read_csv("input/RIN.csv") %>%
    mutate(type = case_when(
      Chain1 == "A" & Chain2 == "A" ~ "ck",
      Chain1 == "B" & Chain2 == "B" ~ "ckr",
      Chain1 == "A" & Chain2 == "B" ~ "inter"
    )) %>% filter(type =="inter") %>%
    filter(class == "xray")
  
  # simplify DF
  rin <- rin %>% select(-PDB, -Chain1, -Chain2, -SS1, -SS2, 
                        -Number.of.atomic.contacts, -class, -type)
  
  # override domain designations for receptor N-terminal Cys and residues
  # coming after
  rin <- rin %>% mutate(dom2 = case_when(
    target_gnccn == "NTr.Cys" ~ "TM1",
    target_gnccn == "NTr.Cp1" ~ "TM1",
    target_gnccn != "NTr.Cp1" & target_gnccn != "NTr.Cp1" ~ dom2
  ))
  
  
  # (2) ADD LABELS
  # add str, unstr labels
  rin <- rin %>% mutate(str_ck = case_when(
    dom1 == "NTc" ~ "NTc",
    dom1 != "NTc" ~ "core_c"
  ))
  
  rin <- rin %>% mutate(str_ckr = case_when(
    dom2 == "NTr" ~ "NTr",
    dom2 == "ECL2" ~ "ECL2",
    dom2 != "NTr" & dom2 != "ECL2" ~ "core_r"
  ))
  # note that ECL1 and ECL3 are short loops and are wrapped into the core,
  # similar to chemokine loops
  
  # (3) GET SUMMARY STATS
  # Get numbers for all complexes individually then average the averages
  
  # all
  all <- rin %>% select(source_gnccn, target_gnccn, str_ck, str_ckr, file) %>% 
    count(str_ck, str_ckr, file)
  all <- all %>% group_by(file) %>% mutate(total_per_file = sum(n)) %>% ungroup()
  all <- all %>% mutate(pct = n / total_per_file)
  all <- all %>% select(-file, -n, -total_per_file) %>% 
    group_by(str_ck, str_ckr) %>% summarise(mean = mean(pct)) %>% ungroup()
  
  # write for cytoscape
  write_csv(all, "output/sse_20190620.csv")
  
