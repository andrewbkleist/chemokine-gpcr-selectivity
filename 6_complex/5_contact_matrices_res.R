# Name:     5_contact_matrices_res.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Supp. Figure 11A,B 

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ggrepel)
  library(bio3d)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_complex/")

##### CONTACT MATRIX ###########################################################
  
  # (1) IMPORT -----------------------------------------------------------------
  
  # Import contacts for all complexes
  
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
  
  
  # (2) RESIDUE ----------------------------------------------------------------
  
  # count occurances of individual contacts
  res <- rin %>% select(source_gnccn, target_gnccn) %>% count(source_gnccn, target_gnccn)
  
  # order positions
  order.ck <- as.factor(c("NTc.Cm10","NTc.Cm9","NTc.Cm8","NTc.Cm7","NTc.Cm6","NTc.Cm5","NTc.Cm4","NTc.Cm3","NTc.Cm2","NTc.Cm1", "NTc.Cm0", "CX.1", "CX.2", "CX.3", "CX.4", "CX.5","cxb1.1","cxb1.2","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.10","cxb1.11","cxb1.14","cxb1.15","cxb1.16","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.4","b1b2.6","b1b2.9","b1b2.10","b1b2.12","b1b2.13","b1b2.14","b1b2.16","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.2","b2b3.4","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10"))
  res$source_gnccn <- factor(res$source_gnccn, levels = rev(order.ck))
  
  order.ckr <- as.factor(c("NTr.Cm8","NTr.Cm7","NTr.Cm6","NTr.Cm5","NTr.Cm4","NTr.Cm3","NTr.Cm2","NTr.Cm1","NTr.Cys","NTr.Cp1","1x24","1x25","1x26","1x27","1x28","1x29","1x30","1x31","1x32","1x33","1x34","1x35","1x36","1x37","1x38","1x39","1x40","1x41","1x42","1x43","1x44","1x45","1x46","1x47","1x48","1x49","1x50","1x51","1x52","1x53","1x54","1x55","1x56","1x57","1x58","1x59","1x60","12x48","12x49","12x50","12x51","2x37","2x38","2x39","2x40","2x41","2x42","2x43","2x44","2x45","2x46","2x47","2x48","2x49","2x50","2x51","2x52","2x53","2x54","2x55","2x56","2x57","2x58","2x59","2x60","2x61","2x62","2x63","2x64","2x65","2x66", "23x50", "3x22","3x23","3x24","3x25","3x26","3x27","3x28","3x29","3x30","3x31","3x32","3x33","3x34","3x35","3x36","3x37","3x38","3x39","3x40","3x41","3x42","3x43","3x44","3x45","3x46","3x47","3x48","3x49","3x50","3x51","3x52","3x53","3x54","3x55","ICL2","4x38","4x39","4x40","4x41","4x42","4x43","4x44","4x45","4x46","4x47","4x48","4x49","4x50","4x51","4x52","4x53","4x54","4x55","4x56","4x57","4x59","4x60","4x61","4x62","4x63", "4x65", "ECL2.Cm8","ECL2.Cm7","ECL2.Cm6","ECL2.Cm5","ECL2.Cm4","ECL2.Cm3","ECL2.Cm2","ECL2.Cm1","45x50","45x51","45x52","ECL2.Cp3","ECL2.Cp4","ECL2.Cp5","ECL2.Cp6","ECL2.Cp7","ECL2.Cp8","ECL2.Cp9", "5x32", "5x34","5x35","5x36","5x37","5x38","5x39","5x40","5x41","5x42","5x43","5x44","5x45","5x46","5x461","5x47","5x48","5x49","5x50","5x51","5x52","5x53","5x54","5x55","5x56","5x57","5x58","5x59","5x60","5x61","5x62","5x63","5x64","5x65","5x66","5x67","5x68","ICL3","6x29","6x30","6x31","6x32","6x33","6x34","6x35","6x36","6x37","6x38","6x39","6x40","6x41","6x42","6x43","6x44","6x45","6x46","6x47","6x48","6x49","6x50","6x51","6x52","6x53","6x54","6x55","6x56","6x57","6x58","6x59","6x60","6x61","6x62","6x63","6x64","ECL3.Cm5","ECL3.Cm4","ECL3.Cm3","ECL3.Cm2","7x23","7x24","7x25","7x26","7x27","7x28","7x29","7x30","7x31","7x32","7x33","7x34","7x35","7x36","7x37","7x38","7x39","7x40","7x41","7x42","7x43","7x45","7x46","7x47","7x48","7x49","7x50","7x51","7x52","7x53","7x54","7x55","7x56","8x47","8x48","8x49","8x50","8x51","8x52","8x53","8x54","8x55","8x56","8x57","8x58","8x59","8x60","8x61","8x62"))
  res$target_gnccn <- factor(res$target_gnccn, levels = order.ckr)
  
  # add descriptive column for number of contacts
  res <- res %>% mutate(no = case_when(
    n == 1 ~ "one",
    n == 2 ~ "two",
    n == 3 ~ "three",
    n == 4 ~ "four"
  ))
  
  # (3) PLOT -------------------------------------------------------------------
  
  # PLOT POSITION MATRIX
  res %>%
    ggplot(aes(target_gnccn, source_gnccn, fill = n)) +
    geom_tile() +
    scale_fill_gradient(low="grey90", high="black") +
    geom_text(aes(target_gnccn, source_gnccn, label = n), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # PLOT HISTO
  res %>%
    ggplot(aes(n)) +
    geom_histogram(binwidth=0.5) +
    theme_minimal()
    
  