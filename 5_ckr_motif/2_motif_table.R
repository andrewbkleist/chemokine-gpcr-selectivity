# Name:     2_motif_analysis.R
# Updated:  20191031
# User:     Andrew Kleist
# Figure:   Figure 5

# Note that this script does not generate any figures but does generate input
# data that contributes to Figure 5

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/5_ckr_motif/")

##### 1:FREQUENCY-OF-OCCURANCE TABLE - NTERM ###################################
  
  # setup data
  data <- read_csv("output/CKR_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
  data <- data %>% select(-motif_no) %>% distinct() 
    # remove doubles of motifs within same sequence
  no.ortho <- read.csv("input/no_ortho.csv")
  no.ortho$protein <- as.character(no.ortho$protein)
  data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$protein)]
  data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$protein)]
  rm(no.ortho)
  data$motif[is.na(data$motif)] <- c("AsnAla") # no gaps are included so all NAs must be AsnAla
  
  # INDIVIDUAL - motif occurance across a single chemokine
  indiv <- data %>% count(motif, protein, class, mer, mask, no_ortho) %>% 
    mutate(pct_ortho = n / no_ortho)
  colnames(indiv)[6] <- c("total_ortho")
  colnames(indiv)[7] <- c("count_ortho")
  indiv <- indiv %>% select(motif, protein, class, mer, mask, count_ortho, total_ortho, pct_ortho)
  
  # FAMILY - motif occurance across each class
    # only considering motifs as occuring in a chemokine if motif appears in
    # ≥50% of orthologous sequences for that chemokine
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
      class == "cc" ~ 10,
      class == "cxc" ~ 6,
      class == "ack" ~ 5,
      class == "xc" ~ 1,
      class == "cx3c" ~ 1
    ))
  freq <- freq %>% mutate(pct_family = count_family / total_family)
  freq$total_super <- 23
  freq <- freq %>% mutate(pct_super = count_super / total_super)
  
  freq <- freq %>% select(motif, protein, class, mer, mask, 
                          count_ortho, total_ortho, pct_ortho,
                          count_family, total_family, pct_family,
                          count_super, total_super, pct_super)
  
  write_csv(freq, "output/CKR_MOTIF_FREQUENCY.csv")
  rm(data, family, freq, indiv, super)
  
  
##### 2:FREQUENCY-OF-OCCURANCE TABLE - ECL2 ####################################
  
  # setup data
  data1 <- read_csv("output/CKR_UNSTRUCTURED_ECL2A_ALL_MERS_CYSLESS.csv")
  data2 <- read_csv("output/CKR_UNSTRUCTURED_ECL2B_ALL_MERS_CYSLESS.csv")
  data <- bind_rows(data1, data2)
  rm(data1, data2)
  data <- data %>% select(-file, -motif_no) %>% distinct() 
  # remove doubles of motifs within same sequence
  
  no.ortho <- read.csv("input/no_ortho.csv")
  no.ortho$protein <- as.character(no.ortho$protein)
  data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$protein)]
  data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$protein)]
  rm(no.ortho)
  data$motif[is.na(data$motif)] <- c("AsnAla") # no gaps are included so all NAs must be AsnAla
  
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
      class == "cc" ~ 10,
      class == "cxc" ~ 6,
      class == "ack" ~ 5,
      class == "xc" ~ 1,
      class == "cx3c" ~ 1
    ))
  freq <- freq %>% mutate(pct_family = count_family / total_family)
  freq$total_super <- 23
  freq <- freq %>% mutate(pct_super = count_super / total_super)
  
  freq <- freq %>% select(motif, protein, class, mer, mask, 
                          count_ortho, total_ortho, pct_ortho,
                          count_family, total_family, pct_family,
                          count_super, total_super, pct_super)
  
  write_csv(freq, "output/CKR_MOTIF_FREQUENCY_ECL2.csv")
  rm(data, family, freq, indiv, super)
  
  
##### 3: MATRIX GENERATION - NTERM #############################################
  
  # import master data
  data <- read_csv("output/CKR_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
  data <- data %>% select(-motif_no) %>% distinct() 
    # remove doubles of motifs within same sequence
  data$motif[is.na(data$motif)] <- c("AsnAla")
  
  # count motif occurances per sequence (all should be 1)
  wide <- data %>% dplyr::count(file, protein, class, motif)
  
  # spread to matrix
  wide <- wide %>% spread(key = motif, value = n)
  wide[is.na(wide)] <- 0
  
  # write csv
  write_csv(wide, "output/CKR_MOTIF_MATRIX.csv")
  
  rm(data, wide)
  
##### 4: MATRIX GENERATION - ECL2 ##############################################
  
  # import master data
  # setup data
  data1 <- read_csv("output/CKR_UNSTRUCTURED_ECL2A_ALL_MERS_CYSLESS.csv")
  data2 <- read_csv("output/CKR_UNSTRUCTURED_ECL2B_ALL_MERS_CYSLESS.csv")
  data <- bind_rows(data1, data2)
  rm(data1, data2)
  data <- data %>% select(-file, -motif_no) %>% distinct() 
  data$motif[is.na(data$motif)] <- c("AsnAla")
  
  # count motif occurances per sequence (all should be 1)
  wide <- data %>% dplyr::count(seqid, protein, class, motif)
  
  # spread to matrix
  wide <- wide %>% spread(key = motif, value = n)
  wide[is.na(wide)] <- 0
  #wide <- wide %>% select(-seqid)
  
  # write csv
  write_csv(wide, "output/CKR_MOTIF_MATRIX_ECL2.csv")
  
  rm(data, wide)
  
  