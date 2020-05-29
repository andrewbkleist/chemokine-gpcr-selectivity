# Name:     0_unstr_str_pie.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Figure 6C 

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_complex/")

  
##### STRUCTURED AND UNSTRUCTURED PIE CHART ####################################

  # (1) IMPORT CONTACTS --------------------------------------------------------
  
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
  
  
  # (2) ADD LABELS -------------------------------------------------------------
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
  
  # (3) GET SUMMARY STATS ------------------------------------------------------
  # Specifically, what is the proportion of chemokine (or receptor) interface
  # positions that are structured or unstructured? This is done in a way that
  # for each chemokine, each position is only counted once. In other words,
  # position CX.1 is only coiunted once for 5UIW even if it makes 3 contacts
  # You want to do this PER STRUCTURE, NORMALIZE BY TOTAL CONTACTS
  
  # chemokine
  ck <- rin %>% select(str_ck, source_gnccn, file) %>% unique() 
  ck.sum <- ck %>% count(file) # get no. interface positions per file
  colnames(ck.sum)[2] <- c("no_pos_per_file")
  ck <- ck %>% count(str_ck, file)
  ck <- left_join(ck, ck.sum)
  ck <- ck %>% mutate(pct = n / no_pos_per_file)
  ck <- ck %>% group_by(str_ck) %>% summarise(mean = mean(pct), sd = sd(pct))
  
  # receptor
  # You want to do this PER STRUCTURE, NORMALIZE BY TOTAL CONTACTS
  
  # chemokine
  ckr <- rin %>% select(str_ckr, target_gnccn, file) %>% unique() 
  ckr.sum <- ckr %>% count(file) # get no. interface positions per file
  colnames(ckr.sum)[2] <- c("no_pos_per_file")
  ckr <- ckr %>% count(str_ckr, file)
  ckr <- left_join(ckr, ckr.sum)
  ckr <- ckr %>% mutate(pct = n / no_pos_per_file)
  ckr <- ckr %>% group_by(str_ckr) %>% summarise(mean = mean(pct), sd = sd(pct))
  
  # (4) PLOT PIE CHART
  ck %>%
    ggplot(aes(x = "", y= mean)) +
    geom_bar(width = 1,size = 1, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  ckr %>%
    ggplot(aes(x = "", y= mean)) +
    geom_bar(width = 1,size = 1, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  