# Name:     3_contact_matrices.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Figure 6F

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_complex/")

##### CONTACT MATRICES #########################################################
  
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
  
  
  # (2) SSE --------------------------------------------------------------------
  
  # You want to calculate number of pairwise domain interactions individually,
  # (ie per complex) THEN get percentages (per complex) to normalize, then
  # get the averages of the normalized pairwise contatcs
  
  # count
  dom <- rin %>% select(dom1, dom2, file) %>% count(dom1, dom2, file)
  dom <- dom %>% group_by(file) %>% mutate(sum = sum(n)) %>% ungroup()
  #dom <- dom %>% mutate(pct = n / sum*100)
  dom <- dom %>% group_by(dom1, dom2) %>% summarize(mean = mean(n), sd = sd(n))
  
  # order domains
  order.ck.dom <- as.factor(unique(c("NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","CX","CX","CX","CX","CX","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","B1","B1","B1","B1","B1","B1","B1","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","B2","B2","B2","B2","B2","B2","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","B3","B3","B3","B3","b3h1","b3h1","b3h1","b3h1","b3h1","b3h1","H","H","H","H","H","H","H","H","H","H")))
  dom$dom1 <- factor(dom$dom1, levels = rev(order.ck.dom))
  order.ckr.dom <- as.factor(unique(c("NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","ICL1","ICL1","ICL1","ICL1","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","ECL1","ECL1","ECL1","ECL1","ECL1","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","ICL3","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","ECL3","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","CT")))
  dom$dom2 <- factor(dom$dom2, levels = order.ckr.dom)
  
  dom <- dom %>% mutate(mean = round(mean, digits = 1))
  
  # PLOT DOMAIN MATRIX
  dom %>%
    ggplot(aes(dom2, dom1, fill = mean)) +
    geom_tile() +
    scale_fill_gradient(low="grey90", high="black") +
    geom_text(aes(dom2, dom1, label = mean), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # write for cyto
  write_csv(dom, "output/dom_cyto.csv")
  