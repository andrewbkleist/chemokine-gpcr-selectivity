# Name:     1_unstr_str_box.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Figure 6D

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_complex/")

##### STRUCTURED AND UNSTRUCTURED BOXPLOTS #####################################
  
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
  
  # (3) COUNT CONTACTS PER ELEMENT PER STRUCTURE
  
  # CHEMOKINE
      # count
      contacts.ck <- rin %>% count(file, ResNum1, str_ck)
      
      # summary stats and p-value
      # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
      ck.sum <- wilcox.test(n ~ str_ck, data = contacts.ck, exact = FALSE)
      ck.sum$p.value
      
      # plot
      contacts.ck %>%
        ggplot(aes(str_ck, n)) +
        geom_boxplot() +
        ylim(0,7) +
        theme_minimal()
      #facet_grid(. ~ file)
  
  # RECEPTOR
    # count
    contacts.ckr <- rin %>% count(file, ResNum2, str_ckr)
    
    # summary stats and p-value
    # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
    ckr.nt_core <- filter(contacts.ckr, str_ckr == "NTr" | str_ckr == "core_r" )
    ckr.nt_core <- wilcox.test(n ~ str_ckr, data = ckr.nt_core, exact = FALSE)
    ckr.nt_core$p.value # significant
    
    ckr.nt_ecl2 <- filter(contacts.ckr, str_ckr == "NTr" | str_ckr == "ECL2" )
    ckr.nt_ecl2 <- wilcox.test(n ~ str_ckr, data = ckr.nt_ecl2, exact = FALSE)
    ckr.nt_ecl2$p.value # not significant
    
    ckr.core_ecl2 <- filter(contacts.ckr, str_ckr == "core_r" | str_ckr == "ECL2" )
    ckr.core_ecl2 <- wilcox.test(n ~ str_ckr, data = ckr.core_ecl2, exact = FALSE)
    ckr.core_ecl2$p.value # significant
    
    # plot
    contacts.ckr %>%
      ggplot(aes(str_ckr, n)) +
      geom_boxplot() +
      ylim(0,7) +
      theme_minimal()
    #facet_grid(. ~ file)
