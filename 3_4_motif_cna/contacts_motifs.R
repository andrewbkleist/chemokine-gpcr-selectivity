# Name:     contacts_motifs.R
# Updated:  20191031
# User:     Andrew Kleist
# Figure:   Figure 3B; Supp. Figure 5A; Figure 5B; Supp. Figure 9A

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(bio3d)
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/3_4_motif_cna/")

##### 1: CK NTERM ##############################################################

  # import contacts, isolate N-termini
  data.res <- read_csv("input/RIN.csv") %>% 
    filter(file =="4rws" | file =="4xt1" | file =="5uiw"| file =="5wb2") %>%
    filter(dom1 == "NTc") %>%
    select(ResNum1, file) %>%
    count(ResNum1, file)

  data.res %>% ggplot(aes(ResNum1, file)) +
    geom_point(aes(size = n), color = "grey40") +
    scale_radius(range = c(1,10)) +
    #scale_color_gradient(low="grey70", high="steelblue4") +
    theme_minimal() +
    coord_flip()
  
##### 2: CKR NTERM #############################################################
  
  # import contacts, isolate N-termini
  data.res <- read_csv("input/RIN.csv") %>% 
    filter(file =="4rws" | file =="4xt1" | file =="5uiw"| file =="5wb2") %>%
    filter(dom2 == "NTr") %>%
    select(ResNum2, file) %>%
    count(ResNum2, file)
  
  data.res %>% ggplot(aes(ResNum2, file)) +
    geom_point(aes(size = n)) +
    scale_radius(range = c(1,10)) +
    #scale_color_gradient(low="grey70", high="steelblue4") +
    theme_minimal() +
    coord_flip()
  
##### 3: CKR ECL2 #############################################################
  
  # import contacts, isolate N-termini
  data.res <- read_csv("input/RIN.csv") %>% 
    filter(file =="4rws" | file =="4xt1" | file =="5uiw"| file =="5wb2") %>%
    filter(dom2 == "ECL2") %>%
    select(target_gnccn, file) %>%
    #filter(target_gnccn != 249) %>%
    count(target_gnccn, file)
  
  order <- c("ECL2.Cm10", "ECL2.Cm9", "ECL2.Cm8", "ECL2.Cm7", "ECL2.Cm6",
             "ECL2.Cm5", "ECL2.Cm4", "ECL2.Cm3", "ECL2.Cm2", "ECL2.Cm1", 
             "45x50", "45x51", "45x52", "ECL2.Cp1", "ECL2.Cp2", "ECL2.Cp3",
             "ECL2.Cp4", "ECL2.Cp5", "ECL2.Cp6", "ECL2.Cp7", "ECL2.Cp8",
             "ECL2.Cp9")
  data.res$target_gnccn <- factor(data.res$target_gnccn, levels = order)
  data.res %>% ggplot(aes(target_gnccn, file)) +
    geom_point(aes(size = n)) +
    scale_radius(range = c(1,10)) +
    #scale_color_gradient(low="grey70", high="steelblue4") +
    theme_minimal() +
    coord_flip()
  
  