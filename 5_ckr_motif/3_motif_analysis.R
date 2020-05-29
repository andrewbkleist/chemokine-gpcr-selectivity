# Name:     3_motif_analysis.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 5C (top),E,F,G,I ; Supp. Figure 9E
  
##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  library(UpSetR)
  library(reshape2)
  library(RColorBrewer)
  library(ggpubr)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/5_ckr_motif/")

##### 0: TIER 1 RANK ORDER #####################################################
# Makes Figure 5C (top)
  
  # import, reformat
  data <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")  
  
  # count the number of appearances of each moif across all 1058 sequences
  data <- data %>% select(motif, count_ortho) %>% group_by(motif) %>% 
    summarize(sum = sum(count_ortho)) %>% ungroup() %>% mutate(pct = sum / 951)
  
  # PLOT
  data$motif <- factor(data$motif, levels = data$motif[order(data$pct)])
  data %>%
    top_n(10, pct) %>%
    ggplot(aes(motif, pct)) +
    geom_point() +
    coord_flip() +
    ylim(0,1) +
    theme_minimal()
  
  # import, reformat
  data <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv")  
  
  # count the number of appearances of each moif across all 1058 sequences
  data <- data %>% select(motif, count_ortho) %>% group_by(motif) %>% 
    summarize(sum = sum(count_ortho)) %>% ungroup() %>% mutate(pct = sum / 951)
  
  # PLOT
  data$motif <- factor(data$motif, levels = data$motif[order(data$pct)])
  data %>%
    top_n(10, pct) %>%
    ggplot(aes(motif, pct)) +
    geom_point() +
    coord_flip() +
    ylim(0,1) +
    theme_minimal()
  
##### 0: HEATMAP - NTERM - BINARIZED ###########################################
# Makes Figure 5E (left)
  
  # import, reformat
  data <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")  

  # select only conserved motifs
  data <- subset(data, data$pct_ortho >= 0.5)

  # select only relevant info (motif, protein, class, pct ortho)
  data <- data %>% select(motif, protein, class, pct_ortho)
  
  # order receptors
  order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                   "ccr5","ccr6","ccr7","ccr8","ccr9",
                                   "ccr10","cxcr1","cxcr2","cxcr3",
                                   "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                   "ackr1","ackr2","ackr3","ackr4",
                                   "ccrl2")))
  
  levels(data$protein)
  data$protein <- factor(data$protein, levels = order.ckr)
  
  # order motifs based on appearance:
  # want most unique to most conserved (global) then most to least conserved among orthologs
  data <- data %>% add_count(motif)
  data <- data %>% arrange(n, protein, -pct_ortho)
  order.motif <- unique(data$motif)
  levels(data$motif)
  data$motif <- factor(data$motif, levels = rev(order.motif))
  
  #1.2 PLOTTING ----------------------------------------------------------------
  
  # MASTER HEATMAP
  data %>%
    ggplot() + 
    geom_raster(aes(protein, motif), fill = "mediumorchid4")+
    theme_minimal() +
    theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  #1.3 COUNT NUMBER OF MOTIFS --------------------------------------------------
  
  test <- data %>% select(motif, n) %>% distinct()
  rm(data, order.ckr, order.motif, test)
  
    
##### 2: HEATMAP - ECL2 - BINARIZED ############################################
# Makes Figure 5E (right)   
  
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv") 
    
    # select only conserved motifs
    data <- subset(data, data$pct_ortho >= 0.5)
    
    # select only relevant info (motif, protein, class, pct ortho)
    data <- data %>% select(motif, protein, class, pct_ortho)
    
    # order receptors
    order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr5","ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    # order motifs based on appearance:
    # want most unique to most conserved (global) then most to least conserved among orthologs
    data <- data %>% add_count(motif)
    data <- data %>% arrange(n, protein, -pct_ortho)
    order.motif <- unique(data$motif)
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.motif))
    
    
    #1.2 PLOTTING ----------------------------------------------------------------
    
    # MASTER HEATMAP
    data %>%
      ggplot() + 
      geom_raster(aes(protein, motif), fill = "mediumorchid4")+
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    #1.3 COUNT NUMBER OF MOTIFS --------------------------------------------------
    
    test <- data %>% select(motif, n) %>% distinct()
    rm(data, order.ckr, order.motif, test)
    
    
##### 3: HEATMAP - ECL2 - PxA AND OTHERS #######################################
# Makes Figure 5I; Supp. Figure 9E  
    
    # (1) ECL2 - PxA
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv")  # 32755 motifs
    
    # select conserved
    data <- data %>% filter(pct_ortho >= 0.5)
    
    # select only relevant info (motif, protein, class, pct ortho)
    data <- data %>% select(motif, protein, class, pct_ortho)
    
    # order receptors
    order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr5","ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    
    # select MOTIF "SLICE"
    data <- data %>% filter(motif == "FP" | motif == "YP" | motif == "FY" | motif == "PY" | motif == "YF" | motif == "PF" | motif == "FxP" | motif == "YxP" | motif == "FxY" | motif == "PxY" | motif == "YxF" | motif == "PxF" |
                              motif == "CK" | motif == "CS" | motif == "AS" | # dummys
                              motif == "AE" | motif == "CG" | motif == "CD" |
                              motif == "CT" | motif == "AD" | motif == "CxT" | 
                              motif == "DV" | motif == "DK" | motif == "KV" |
                              motif == "GL" | motif == "DF" | motif == "EQ")
    
    order.mot <- as.factor(c("FP", "FY", "PF", "PY", "YF", "YP", "FxP", "FxY", "PxF", "PxY", "YxF", "YxP"))
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.mot))
    
    data %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumorchid4") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    

    # (2) SULFOTYROSINE MOTIFS
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")  # 32755 motifs
    
    # select conserved
    data <- data %>% filter(pct_ortho >= 0.5)
    
    # select only relevant info (motif, protein, class, pct ortho)
    data <- data %>% select(motif, protein, class, pct_ortho)
    
    # order receptors
    order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr5","ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    # select MOTIF "SLICE"
    data <- data %>% filter( motif == "DxT" |
                              motif == "FP" | motif == "YP" | motif == "FY" | # dummys
                               motif == "PY" | motif == "YF" | motif == "PF" | 
                               motif == "FxP" | motif == "YxP" | motif == "FxY" | 
                               motif == "PxY" | motif == "YxF" | motif == "PxF" |
                              motif == "CK" | motif == "CS" | motif == "AS" | 
                              motif == "AE" | motif == "CG" | motif == "PC" |
                              motif == "CT" | motif == "DY")
    
    order.mot <- as.factor(c("DxT", "DY", "PC", "PF", "PY", "YF", "YP", "FxP", "FxY", "PxF", "PxY", "YxF", "YxP"))
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.mot))
    
    data %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumorchid4") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    
##### 3: CONSERVATION VS. UNIQUENESS BUBBLE PLOT ###############################
# Makes Figure 5F,G
    
    # import, reformat
    nt <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")  
    ecl2 <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv")  
    
    # select only conserved motifs
    nt <- nt %>% filter(pct_ortho >= 0.5)
    ecl2 <- ecl2 %>% filter(pct_ortho >= 0.5)
    
    # add column for uniqueness
    nt <- nt %>% mutate(uniqueness = 1- pct_super + (1/23)) # uniqueness = 1-[5/23] + (1/23)
    ecl2 <- ecl2 %>% mutate(uniqueness = 1- pct_super + (1/23)) # uniqueness = 1-[5/23] + (1/23)
    
    # RECLASSIFY PEPTIDE LENGTH CORRECTING FOR GAPS
    nt <- nt %>% mutate(meradj = case_when(
      mer == "mer2" ~ "2",
      grepl('x', motif) & mer == "mer3" ~ "2",
      !grepl('x', motif) & mer == "mer3" ~ "3",
      !grepl('x', motif) & mer == "mer3" ~ "3",
      !grepl('x', motif) & mer == "mer4" ~ "4",
      grepl("[A-Z]x[A-Z][A-Z]", motif)  ~ "3",
      grepl("[A-Z][A-Z]x[A-Z]", motif)  ~ "3",
      grepl("[A-Z]xx[A-Z]", motif)  ~ "2"
    ))
    ecl2 <- ecl2 %>% mutate(meradj = case_when(
      mer == "mer2" ~ "2",
      grepl('x', motif) & mer == "mer3" ~ "2",
      !grepl('x', motif) & mer == "mer3" ~ "3",
      !grepl('x', motif) & mer == "mer3" ~ "3",
      !grepl('x', motif) & mer == "mer4" ~ "4",
      grepl("[A-Z]x[A-Z][A-Z]", motif)  ~ "3",
      grepl("[A-Z][A-Z]x[A-Z]", motif)  ~ "3",
      grepl("[A-Z]xx[A-Z]", motif)  ~ "2"
    ))
    ecl2$sse <- c("ecl2") 
    nt$sse <- c("nt")
    master <- bind_rows(ecl2, nt)
    
    # PLOT ALL (15561 motifs)
    nt %>%
      ggplot(aes(uniqueness, pct_ortho, size=meradj, alpha=0.9)) +
      geom_jitter(color = "mediumorchid4") +
      xlim(0,1) +
      ylim(0.49,1) +
      theme_minimal()
   
    ecl2 %>%
      ggplot(aes(uniqueness, pct_ortho, size=meradj, alpha=0.9)) +
      geom_jitter(color = "mediumorchid4") +
      xlim(0,1) +
      ylim(0.49,1) +
      theme_minimal()

    # specific motif
    nt %>%
      filter(motif == "DY") %>%
      filter(protein == "ccr5") %>%
      ggplot(aes(uniqueness, pct_ortho, size=meradj)) +
      geom_point() +
      xlim(0,1) +
      ylim(0.49,1) +
      #scale_color_brewer(palette="Blues") +
      theme_minimal()
    
    ecl2 %>%
      filter(motif == "HFP") %>%
      filter(protein == "ccr5") %>%
      ggplot(aes(uniqueness, pct_ortho, size=meradj)) +
      geom_point() +
      xlim(0,1) +
      ylim(0.49,1) +
      #scale_color_brewer(palette="Blues") +
      theme_minimal()
    
    
    rin <- read_csv("../3A/input/RIN.csv")
    
    
    # NOW PLOT COMPARISON OF UNIQUENESS
    order.mot <- as.factor(c("nt", "ecl2"))
    master$sse <- factor(master$sse, levels = order.mot)
    
    master %>%
      filter(count_super >1) %>%
      ggplot(aes(uniqueness, fill = sse, color = sse)) +
      geom_histogram(position="dodge", binwidth = .05) +
      geom_freqpoly(binwidth = .05) +
      scale_fill_manual(values = c("grey70", "grey40")) +
      scale_color_manual(values = c("grey70", "grey40")) +
      #geom_density() +
      #xlim(0.3,1) +
      geom_vline( xintercept = 0.9565217) + # shared in 2
      geom_vline( xintercept = 0.826087) + # shared in 5
      geom_vline( xintercept = 0.6086957) + # shared in 10
      theme_minimal()
    
    # summary stats and p-value
    # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
    test <- wilcox.test(uniqueness ~ sse, data = master, exact = FALSE)
    test$p.value
    
##### 4: RECEPTOR-SPECIFIC MOTIF MATRICES #####################################
# Makes Supp. Figure 10
    
    # (4.1) CXCR4 N-TERM - 2-MER -----------------------------------------------
    
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")
    data <- data %>% filter(pct_ortho >= 0.5) 
    
    
    # MASTER HEATMAP - MOTIF "SLICE"
    data <- data %>% filter(motif == "ME" | motif == "EG" | motif == "GI" | 
                              motif == "IS" | motif == "SI" | motif == "IY" | 
                              motif == "YT" | motif == "TS" | motif == "SD" | 
                              motif == "DN" | motif == "NY" | motif == "YT" |
                              motif == "TE" | motif == "EE" | motif == "EM" |
                              motif == "MG" | motif == "GS" | motif == "SG" |
                              motif == "GD" | motif == "DY" | motif == "YD" |
                              motif == "DS" | motif == "FS" | motif == "SM" |
                              motif == "MK" | motif == "KE" | motif == "EP" |
                              motif == "PC" | motif == "CF")
    
    order.mot <- as.factor(unique(c("ME","EG","GI","IS","SI","IY","YT","TS","SD","DN","NY","YT","TE","EE","EM","MG","GS","SG","GD","DY","YD","DS","FS","SM","MK","KE","EP","PC","CF")))
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.mot))
    
    # order chemokines
    order.ckr <- as.factor(tolower(c("cxcr4", "ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr5","ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    data %>%
      #filter(pct_ortho > 0.5) %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumorchid4") +
      #scale_fill_gradient(low = "white", high = "steelblue3") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    # (4.2) CXCR4 ECL2 - 2-MER -------------------------------------------------
    
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv")  
    data <- data %>% filter(pct_ortho >= 0.5) 
    
    # MASTER HEATMAP - MOTIF "SLICE"
    data <- data %>% filter(motif == "NV" | motif == "VS" | motif == "SE" | 
                              motif == "EA" | motif == "AD" | motif == "DD" | 
                              motif == "DR" | motif == "RY" | motif == "YI" | 
                              motif == "IC" | motif == "CD" | motif == "DR" |
                              motif == "RF" | motif == "FY" | motif == "YP" |
                              motif == "PN" |
                              motif == "CS" | motif == "CG" | motif == "CK" | #dummy motifs
                              motif == "CE" | motif == "CR" | motif == "CQ" |
                              motif == "CS" | motif == "FP" | motif == "FL" | 
                              motif == "CT" | motif == "DS" | motif == "ER" |
                              motif == "DV" | motif == "DK" | motif == "GL" |
                              motif == "DF" )
    
    order.mot <- as.factor(unique(c("NV","VS","SE","EA","AD","DD","DR","RY","YI","IC","CD","DR","RF","FY","YP","PN")))
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.mot))
    
    # order chemokines
    order.ckr <- as.factor(tolower(c("cxcr4","ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr5","ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    data %>%
      #filter(pct_ortho > 0.5) %>%
      ggplot() + 
      geom_tile(aes(protein, motif),  fill = "mediumorchid4") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    # (4.3) CCR5 N-TERM - 2-MER ------------------------------------------------
    
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")
    data <- data %>% filter(pct_ortho >= 0.5) 
    
    
    # MASTER HEATMAP - MOTIF "SLICE"
    data <- data %>% filter(motif == "MD" | motif == "DY" | motif == "YQ" | 
                              motif == "QV" | motif == "VS" | motif == "SS" | 
                              motif == "SP" | motif == "PI" | motif == "IY" | 
                              motif == "YD" | motif == "DI" | motif == "IN" |
                              motif == "NY" | motif == "YY" | motif == "YT" |
                              motif == "TS" | motif == "SE" | motif == "EP" |
                              motif == "PC" | motif == "CQ")
    
    order.mot <- as.factor(unique(c("MD", "DY", "YQ", "QV", "VS", "SS", "SP", "PI", "IY", "YD", "DI", "IN", "NY", "YY", "YT", "TS", "SE", "EP", "PC", "CQ")))
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.mot))
    
    # order chemokines
    order.ckr <- as.factor(tolower(c("ccr5", "ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    data %>%
      #filter(pct_ortho > 0.5) %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumorchid4") +
      #scale_fill_gradient(low = "white", high = "steelblue3") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    # (4.4) CCR5 ECL2 - 2-MER --------------------------------------------------
    
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY_ECL2.csv")  
    data <- data %>% filter(pct_ortho >= 0.5) 
    
    # MASTER HEATMAP - MOTIF "SLICE"
    data <- data %>% filter(motif == "RS" | motif == "SQ" | motif == "QK" | 
                              motif == "KE" | motif == "EG" | motif == "GL" | 
                              motif == "LH" | motif == "HY" | motif == "YT" | 
                              motif == "TC" | motif == "CS" | motif == "SS" |
                              motif == "SH" | motif == "HF" | motif == "FP" |
                              motif == "PY" |
                              motif == "YS" |
                            motif == "NV" | motif == "VS" | motif == "SE" | 
                              motif == "EA" | motif == "AD" | motif == "DD" | 
                              motif == "DR" | motif == "RY" | motif == "YI" | 
                              motif == "IC" | motif == "CD" | motif == "DR" |
                              motif == "RF" | motif == "FY" | motif == "YP" |
                              motif == "PN" |
                              motif == "CS" | motif == "CG" | motif == "CK" | #dummy motifs
                              motif == "CE" | motif == "CR" | motif == "CQ" |
                              motif == "CS" | motif == "FP" | motif == "FL" | 
                              motif == "CT" | motif == "DS" | motif == "ER" |
                              motif == "DV" | motif == "DK" | motif == "GL" |
                              motif == "DF")
    
    order.mot <- as.factor(unique(c("RS", "SQ", "QK", "KE","EG", "GL", "LH", "HY", "YT", "TC", "CS", "SS", "SH" ,"HF", "FP", "PY", "YS", "CG", "CK", "CE", "CR", "CQ", "CS", "FP", "FL", "CT")))
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.mot))
    
    # order chemokines
    order.ckr <- as.factor(tolower(c("ccr5", "ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr4", "cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    data %>%
      #filter(pct_ortho > 0.5) %>%
      ggplot() + 
      geom_tile(aes(protein, motif),  fill = "mediumorchid4") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    
    # (4.5) SULFOTYROSINE ------------------------------------------------------
    
    # import, reformat
    data <- read_csv("output/CKR_MOTIF_FREQUENCY.csv")  
    data <- data %>% filter(pct_ortho >= 0.5) 

    data <- data %>% filter(motif == "DY" | motif == "YD" | 
                              motif == "EY" | motif == "YE" | 
                              motif == "DxY" | motif == "DxxY" | 
                              motif == "YxxD" | motif == "YxD" |
                              motif == "ExY" | motif == "ExxY" | 
                              motif == "YxxE" | motif == "YxE")

    order.mot <- as.factor(unique(c("DY","YD", "EY", "YE",
                                    "DxY","YxD", "ExY", "YxE",
                                    "DxxY", "YxxD", "ExxY", "YxxE")))
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.mot))

    # order chemokines
    order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                     "ccr5","ccr6","ccr7","ccr8","ccr9",
                                     "ccr10","cxcr1","cxcr2","cxcr3",
                                     "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                     "ackr1","ackr2","ackr3","ackr4",
                                     "ccrl2")))
    
    levels(gluasp$protein)
    gluasp$protein <- factor(gluasp$protein, levels = order.ckr)
    data$protein <- factor(data$protein, levels = order.ckr)
    
    data %>%
      filter(pct_ortho >= 0.5) %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumorchid4") +
      #scale_fill_gradient(low = "white", high = "steelblue3") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

    