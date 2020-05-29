# Name:     3_get_motifs.R
# Updated:  20191026
# Author:   Andrew Kleist
# Figure:   Figure 6 (makes input data)

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################
  
  library(tidyverse)
  library(ggrepel)
  library(ggforce)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")

##### 1: CHEMOKINE (20191024) ##################################################
  
  # import MOTIFS
  lookup.motif <- read_csv("input/20191218_motif_conversion_3mer.csv")
  
  # import MOTIF CONSERVATION
  motif.cons <- read_csv("../3_ck_motif/output/CK_MOTIF_FREQUENCY.csv") %>% 
    #filter(pct_ortho >= 0.5) %>% 
    filter(mer == "mer3")
  motif.cons <- motif.cons %>% dplyr::select(motif, protein, pct_ortho)
  motif.cons$protein <- toupper(motif.cons$protein)
  
  # JOIN ...
  lookup.motif <- left_join(lookup.motif, motif.cons)
  lookup.motif <- lookup.motif %>% filter(!is.na(motif))
  lookup.motif <- lookup.motif %>% mutate(tier = case_when(
    pct_ortho >= 0.5 ~ "motif",
    pct_ortho < 0.5 ~ "fragment"
  ))
  rm(motif.cons)
  write_csv(lookup.motif, "output/ck_motif_conservation.csv") 
  

##### 2: RECEPTOR - NTERM (20191026) ###########################################
  
  # import MOTIFS
  lookup.motif.nt <- read_csv("input/20191218_motif_conversion_3mer_ckr_nterm.csv")
  
  # import MOTIF CONSERVATION
  motif.cons.nt <- read_csv("../5_ckr_motif/output/CKR_MOTIF_FREQUENCY.csv") %>% 
    #filter(pct_ortho >= 0.5) %>% 
    filter(mer == "mer3")
  motif.cons.nt <- motif.cons.nt %>% dplyr::select(motif, protein, pct_ortho)
  motif.cons.nt$protein <- toupper(motif.cons.nt$protein)
  
  # JOIN ...
  lookup.motif.nt <- left_join(lookup.motif.nt, motif.cons.nt)
  lookup.motif.nt <- lookup.motif.nt %>% filter(!is.na(motif))
  lookup.motif.nt <- lookup.motif.nt %>% mutate(tier = case_when(
    pct_ortho >= 0.5 ~ "motif",
    pct_ortho < 0.5 ~ "fragment"
  ))
  rm(motif.cons.nt)
  write_csv(lookup.motif.nt, "output/ckr_nt_motif_conservation.csv")  

  
##### 3: RECEPTOR - ECL2 (20191026) ############################################
  
  # import MOTIFS
  lookup.motif.ecl2 <- read_csv("input/20191218_motif_conversion_3mer_ckr_ecl2.csv")

  # import MOTIF CONSERVATION
  # recall that ECL2a and ECL2b are combined into one data frame
  motif.cons.ecl2 <- read_csv("../5_ckr_motif/output/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
    #filter(pct_ortho >= 0.5) %>% 
    filter(mer == "mer3")
  motif.cons.ecl2 <- motif.cons.ecl2 %>% dplyr::select(motif, protein, pct_ortho)
  motif.cons.ecl2$protein <- toupper(motif.cons.ecl2$protein)
  
  # JOIN ...
  lookup.motif.ecl2 <- left_join(lookup.motif.ecl2, motif.cons.ecl2)
  lookup.motif.ecl2 <- lookup.motif.ecl2 %>% filter(!is.na(motif))
  lookup.motif.ecl2 <- lookup.motif.ecl2 %>% mutate(tier = case_when(
    pct_ortho >= 0.5 ~ "motif",
    pct_ortho < 0.5 ~ "fragment"
  ))
  rm(motif.cons.ecl2)
  write_csv(lookup.motif.ecl2, "output/ckr_ecl2_motif_conservation.csv")  
  