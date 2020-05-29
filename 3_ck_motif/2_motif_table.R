# Name:     2_motif_analysis.R
# Updated:  20191031
# User:     Andrew Kleist
# Figure:   Figure 3

# Note that this script does not generate any figures but does generate input
# data that contributes to Figure 3

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/3_ck_motif/")

##### 1:FREQUENCY-OF-OCCURANCE TABLE ###########################################
  
  # (1) FILTERING & CLEANING
  data <- read_csv("output/CK_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
  
  # REMOVE DOUBLES
  data <- data %>% select(-motif_no) %>% distinct() 
  
  # ADD NUMBER OF ORTHOLOGS
  no.ortho <- read.csv("input/no_ortho.csv")
  data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$ck)]
  data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$ck)]
  rm(no.ortho)
  
  # REPLACE "NA" WITH ASN-ALA
  data$motif[is.na(data$motif)] <- c("AsnAla")
  
  # REMOVE XCL1, XCL2, CCL4L1
  data <- data %>% filter(protein != "xcl1" & protein != "xcl2" & protein != "ccl4l1")
  
  # (2) CONSERVATION SCORING AT DIFFERENT LEVELS (ORTHOLOGS, FAMILY, ETC)
  # INDIVIDUAL - motif occurance across a single chemokine
  indiv <- data %>% count(motif, protein, class, mer, mask, no_ortho) %>% 
    mutate(pct_ortho = n / no_ortho)
  colnames(indiv)[6] <- c("total_ortho")
  colnames(indiv)[7] <- c("count_ortho")
  indiv <- indiv %>% select(motif, protein, class, mer, mask, count_ortho, total_ortho, pct_ortho)
  
  # FAMILY - motif occurance across each class
    # only considering motifs as occuring in a chemokine if motif appears in
    # >50% of orthologous sequences for that chemokine
    # Evaluation occurs on paralog level, not abundance of ortholog sequences
    # due to assymetric data sets for each chemokine (ie some have 50 sequences)
  family <- indiv %>% filter(pct_ortho >= 0.5) %>% count(motif, class)
  colnames(family)[3] <- c("count_family")
    
  freq <- left_join(indiv, family)
  freq[is.na(freq)] <- 0
  
  # SUPERFAMILY - motif occurance across all chemokines
  super <- indiv %>% filter(pct_ortho >= 0.5) %>% count(motif)
  colnames(super)[2] <- c("count_super")
  
  freq <- left_join(freq, super)
  freq[is.na(freq)] <- 0
  
  # ADD PCT COLUMNS
  freq <- freq %>%
    mutate(total_family = case_when(
      class == "cc" ~ 27, # recall that CCL4L1 has been removed, so 27 instead of 28
      class == "cxc" ~ 17,
      class == "cx3c" ~ 1 # note that XC family has been removed
    ))
  
  freq <- freq %>% mutate(pct_family = count_family / total_family)
  freq$total_super <- 43
  freq <- freq %>% mutate(pct_super = count_super / total_super)
  
  # (3) CLEAN AND WRITE
  freq <- freq %>% select(motif, protein, class, mer, mask, 
                          count_ortho, total_ortho, pct_ortho,
                          count_family, total_family, pct_family,
                          count_super, total_super, pct_super)
  
  write_csv(freq, "output/CK_MOTIF_FREQUENCY.csv")
  rm(data, family, freq, indiv, super)
  
##### 2: MATRIX GENERATION #####################################################
  
  # FILTERING & CLEANING
  data <- read_csv("output/CK_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
  
  # REMOVE DOUBLES
  data <- data %>% select(-motif_no) %>% distinct() 
  
  # REPLACE "NA" WITH ASN-ALA
  data$motif[is.na(data$motif)] <- c("AsnAla")
  
  # REMOVE XCL1, XCL2, CCL4L1
  data <- data %>% filter(protein != "xcl1" & protein != "xcl2" & protein != "ccl4l1")
  
  # count motif occurances per sequence (all should be 1)
  wide <- data %>% dplyr::count(file, protein, class, motif)
  
  # spread to matrix
  wide <- wide %>% spread(key = motif, value = n)
  wide[is.na(wide)] <- 0
  
  # write csv
  write_csv(wide, "output/CK_MOTIF_MATRIX.csv")
  
  rm(data, wide)
  
  