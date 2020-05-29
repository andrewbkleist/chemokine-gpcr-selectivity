# Name:     3_motif_analysis.R
# Updated:  20191031
# User:     Andrew Kleist
# Figures:  Figure 3F,G, H; Supp. Figure 5B (left) 

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  library(UpSetR)
  library(reshape2)
  library(RColorBrewer)
  library(ggpubr)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/3_ck_motif/")

##### 0: TIER 1 RANK ORDER (updt: 20190618) ####################################
  
  # import, reformat
  data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")  
  
  # count the number of appearances of each moif across all 1058 sequences
  data <- data %>% select(motif, count_ortho) %>% group_by(motif) %>% 
    summarize(sum = sum(count_ortho)) %>% ungroup() %>% mutate(pct = sum / 1058)
  
  # PLOT
  data$motif <- factor(data$motif, levels = data$motif[order(data$pct)])
  data %>%
    top_n(20, pct) %>%
    ggplot(aes(motif, pct)) +
    geom_point() +
    coord_flip() +
    ylim(0,1) +
    theme_minimal()
  
##### 1: HEATMAP - BINARIZED (updt: 20190523) ##################################
# Makes Figure 3F
  
  # 1.1 SETUP ------------------------------------------------------------------
    # For every motif and every chemokine, give the percent ortholog score of 
    # that motif, colored by class and intensity
    
    # import, reformat
    data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")  
  
    # count numbers in each quad
    count.quad <- data %>% mutate(quad = case_when(
      pct_ortho >= 0.5 & count_super == 1 ~ "TR",
      pct_ortho >= 0.5 & count_super > 1 ~ "TL",
      pct_ortho < 0.5 ~ "B"
    ))
    
    count.quad <- count.quad %>% count(quad)
    
    # select only conserved motifs
    data <- data %>% filter(pct_ortho >= 0.5)
  
    # select only relevant info (motif, protein, class, pct ortho)
    data <- data %>% select(motif, protein, class, pct_ortho)
    
    # order chemokines
    order.ck <- as.factor(tolower(c("ccl1","ccl2", "ccl3", "ccl3l1",
                                    "ccl4","ccl4l1","ccl5","ccl7","ccl8",
                                    "ccl11","ccl13","ccl14","ccl15",
                                    "ccl16","ccl17","ccl18","ccl19",
                                    "ccl20","ccl21","ccl22","ccl23",
                                    "ccl24","ccl25","ccl26","ccl27",
                                    "ccl28","cxcl1","cxcl2","cxcl3",
                                    "cxcl4","cxcl4l1","cxcl5","cxcl6",
                                    "cxcl7","cxcl8","cxcl9","cxcl10",
                                    "cxcl11","cxcl12","cxcl13", "cxcl14", "cxcl16",
                                    "cxcl17", "cx3cl1", "xcl1", "xcl2")))
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ck)
    
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
      geom_tile(aes(protein, motif), fill = "mediumorchid4")+
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    # ELR HEATMAP (adding additional motifs to get full complement of chemokines on x-axis)
    data %>%
      filter(motif == "ELR") %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumpurple4")+
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
##### 3: CONSERVATION VS. UNIQUENESS BUBBLE PLOT ###############################
# Makes Figure 3H
    
    # import, reformat
    data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")  

    # select only conserved motifs
    data <- data %>% filter(pct_ortho >= 0.5)
    
    # add uniqueness column
    data <- data %>% mutate(uniqueness = 1- pct_super + (1/43)) # uniqueness = 1 - [pct_super] + (1/43)
    
    # count numbers of unique and shared
    how.many <- data %>% select(motif,count_super) %>% distinct()
    
    # RECLASSIFY PEPTIDE LENGTH CORRECTING FOR GAPS
    data <- data %>% mutate(meradj = case_when(
      mer == "mer2" ~ "2",
      grepl('x', motif) & mer == "mer3" ~ "2",
      !grepl('x', motif) & mer == "mer3" ~ "3",
      !grepl('x', motif) & mer == "mer3" ~ "3",
      !grepl('x', motif) & mer == "mer4" ~ "4",
      grepl("[A-Z]x[A-Z][A-Z]", motif)  ~ "3",
      grepl("[A-Z][A-Z]x[A-Z]", motif)  ~ "3",
      grepl("[A-Z]xx[A-Z]", motif)  ~ "2"
    ))
    
    # CONSERVATION AND UNIQUENESS BUBBLE PLOT - PLOT ALL (15561 motifs)
    data %>%
      ggplot(aes(uniqueness, pct_ortho, size=meradj, alpha=0.9)) +
      geom_jitter(color = "mediumorchid4") +
      xlim(0.79,1) +
      geom_vline(xintercept = 0.9767442) +
      # ylim(0,1) +
      #scale_color_brewer(palette="Blues") +
      theme_minimal()
  
    # BUBBLE PLOT SLICE HISTOGRAM - UNIQUENESS
    data %>%
      ggplot(aes(uniqueness)) +
      geom_histogram() +
      #ylim(0,2500) +
      theme_minimal()
    
    # BUBBLE PLOT SLICE HISTOGRAM - CONSERVATION
    data %>%
      ggplot(aes(pct_ortho)) +
      geom_histogram() +
      #ylim(0,2500) +
      theme_minimal() +
      facet_grid(meradj ~ .)
    
    # CONSERVATION AND UNIQUENESS BUBBLE PLOT - PLOT ONE
    data %>%
      filter(motif == "ILP") %>%
      filter(protein == "ccl28") %>%
      ggplot(aes(uniqueness, pct_ortho, size=meradj, alpha=0.9, color=meradj)) +
      geom_point() +
      xlim(0.8,1) +
      ylim(0.5,1) +
      scale_color_brewer(palette="Blues") +
      theme_minimal()
    
##### 4: FUNCTION VS. LENGTH (UPDATED 20190524) ################################
# Makes Figure 3G
    
    # (3.1) ADD CYS-INDEX ------------------------------------------------------
    
      # import, reformat
      data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")
    
      # add column for uniqueness
      data <- data %>% mutate(uniqueness = 1- pct_super + (1/43)) # uniqueness = 1 - [pct_super] + (1/43)
    
      # examine only CONSERVED motifs
      data <- data %>% filter(pct_ortho >= 0.5)
      
      # RECLASSIFY PEPTIDE LENGTH CORRECTING FOR GAPS
      data <- data %>% mutate(meradj = case_when(
        mer == "mer2" ~ "2",
        grepl('x', motif) & mer == "mer3" ~ "2",
        !grepl('x', motif) & mer == "mer3" ~ "3",
        !grepl('x', motif) & mer == "mer3" ~ "3",
        !grepl('x', motif) & mer == "mer4" ~ "4",
        grepl("[A-Z]x[A-Z][A-Z]", motif)  ~ "3",
        grepl("[A-Z][A-Z]x[A-Z]", motif)  ~ "3",
        grepl("[A-Z]xx[A-Z]", motif)  ~ "2"
      ))
      
      # GET HUMAN CYSTEINE INDEX
      # Note that by aggregating motifs to calculate conservation
      # across orthologs, you lose the index information from each 
      # individual ortholog - to re-map, you will have to thus
      # choose one ortholog (human) to map to, even if this does
      # not represent the original index of that particular ortholog,
      # which may have been shifted
      
      # (1) import, etc
      raw <- read_csv("output/CK_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
      raw.hum <- subset(raw, grepl('HUMAN', file) | grepl('human', file))
      raw.hum$motif[is.na(raw.hum$motif)] <- c("AsnAla")
      
      # (2) add indices to "data" df
      index <- raw.hum %>% select(protein, motif, motif_no)
      data <- left_join(data, index)
        # note that 166 motifs are NAs, likely because they are conserved motifs
        # that do not appear in humans
      data <- data %>% filter(!is.na(motif_no))
      rm(index)
      
      # (3) Add motif distance from CX.1
      # Note that current indices are from the "front" of the sequence
      # but that you need to adjust to the "back" of the sequence; all of these
      # will vary based upon the length of the peptide
      mer <- raw.hum %>% group_by(protein, mer) %>% 
        filter(motif_no == max(motif_no)) %>%
        select(protein, mer, motif_no) %>%
        distinct()
      colnames(mer)[3] <- c("max_motif")
      data <- left_join(data, mer)
      
      # (4) calculate maximum motif
      data <- data %>% mutate(motif_adj = case_when(
        mer == "mer2" ~ motif_no - max_motif - 2,
        mer == "mer3" ~ motif_no - max_motif - 3,
        mer == "mer4" ~ motif_no - max_motif - 4
        ))
      
      # (5) annotate
      # make character
      data <- data %>% mutate(cys_adj_pos_chr = case_when(
        motif_adj == -12 ~ "m12",
        motif_adj == -11 ~ "m11",
        motif_adj == -10 ~ "m10",
        motif_adj == -9 ~ "m9",
        motif_adj == -8 ~ "m8",
        motif_adj == -7 ~ "m7",
        motif_adj == -6 ~ "m6",
        motif_adj == -5 ~ "m5",
        motif_adj == -4 ~ "m4",
        motif_adj == -3 ~ "m3",
        motif_adj == -2 ~ "m2",
        motif_adj == -1 ~ "m1"      
      ))
      
      # note we are now going to exclude positions that are farther away than -12
      # (average N-terminal length = 12 residues)
      
      # (6) plot
      order <- c("m12", "m11", "m10", "m9", "m8", "m7", "m6", "m5", "m4",
                 "m3", "m2", "m1")
      levels(data$cys_adj_pos_chr)
      data$cys_adj_pos_chr <- factor(data$cys_adj_pos_chr, levels = order)
      
      # CONSERVATION BY LENGTH
      data %>%
        #filter(mer == "mer2") %>%
        filter(!is.na(cys_adj_pos_chr)) %>%
        ggplot() +
        #geom_line() +
        geom_smooth(aes(motif_adj, pct_ortho)) +
        #geom_jitter(aes(cys_adj_pos, uniqueness)) +
        scale_x_continuous( limits = c(-12 , -2), breaks = c(-12:-2)) +
        #xlim(-12 , -2) +
        coord_flip() +
        theme_minimal()
      
      # UNIQUENESS BY LENGTH
      #http://r-statistics.co/Loess-Regression-With-R.html
      data %>%
        filter(!is.na(cys_adj_pos_chr)) %>%
        #filter(mer == "mer2") %>%
        ggplot() +
        #geom_line() +
        geom_smooth(aes(motif_adj, uniqueness)) +
        #geom_jitter(aes(cys_adj_pos, uniqueness)) +
        scale_x_continuous(limits = c(-12 , -2), breaks = c(-12:-2)) +
        #xlim(-12,-2) +
        coord_flip() +
        theme_minimal() 
      
      # SIGNFICANCE TESTING
      pairwise.wilcox.test(data$uniqueness, data$cys_adj_pos_chr,
                           p.adjust.method = "BH")
      pairwise.wilcox.test(data$count_super, data$cys_adj_pos_chr,
                           p.adjust.method = "BH")
      pairwise.wilcox.test(data$pct_ortho, data$cys_adj_pos_chr,
                           p.adjust.method = "BH")
            
##### 5: CHEMOKINE-SPECIFIC MOTIF MATRICES #####################################
      
      # (5.1) CCL5 - 2-MER -----------------------------------------------------
      
      # import, reformat
      data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")  
      data <- data %>% filter(pct_ortho >= 0.5)
      
      # # convert to matrix, fill NA with zero, convert back to long
      # data <- data %>% select(motif, protein, pct_ortho) %>% spread(protein, pct_ortho)
      # data[is.na(data)] <- 0
      # data <- data %>% gather(protein, pct_ortho, 2:ncol(data))
      
      # MASTER HEATMAP - MOTIF "SLICE"
      data <- data %>% filter(motif == "SP" | motif == "PY" | motif == "YS" | 
                                motif == "SS" | motif == "SD" | motif == "DT" | 
                                motif == "TT" | motif == "TP" |
                                motif == "PxA" | motif == "PxS" | motif == "PxT" |
                                motif == "AA" | motif == "AS" | motif == "SA" |    # "dummy motifs so that all chemokines are shown
                                motif == "AT" | motif == "EL" | motif == "LR" |    # "dummy motifs so that all chemokines are shown
                                motif == "LS" | motif == "LL" | motif == "VA" |    # "dummy motifs so that all chemokines are shown
                                motif == "SPY" | motif == "PYS" | motif == "YSS" | # "dummy motifs so that all chemokines are shown
                                motif == "SSD" | motif == "SDT" | motif == "DTT" | # "dummy motifs so that all chemokines are shown
                                motif == "TTP" |                                   # "dummy motifs so that all chemokines are shown
                                motif == "SP" | motif == "PY" | motif == "YS" |    # "dummy motifs so that all chemokines are shown
                                motif == "SS" | motif == "SD" | motif == "DT" |    # "dummy motifs so that all chemokines are shown
                                motif == "TT" | motif == "TP" |                    # "dummy motifs so that all chemokines are shown
                                motif == "KPV" | motif == "PVS" | motif == "VSL" | # "dummy motifs so that all chemokines are shown
                                motif == "SLS" | motif == "LSY" | motif == "SYR" | # "dummy motifs so that all chemokines are shown
                                motif == "AR" | motif == "AE" |motif == "SLS" |    # "dummy motifs so that all chemokines are shown
                                motif == "ED" |motif == "SLS" | motif == "AK" |    # "dummy motifs so that all chemokines are shown
                                motif == "GR" | motif == "FK" | motif == "LE" |    # "dummy motifs so that all chemokines are shown
                                motif == "SK" | motif == "EG" | motif == "GV" )    # "dummy motifs so that all chemokines are shown
      
      order.mot <- as.factor(c("SP","PY","YS","SS","SD","DT","TT", "TP"))
      levels(data$motif)
      data$motif <- factor(data$motif, levels = order.mot)
      
      # order chemokines
      order.ck <- as.factor(tolower(c("ccl5", "ccl1","ccl2", "ccl3", "ccl3l1",
                                      "ccl4","ccl4l1","ccl7","ccl8",
                                      "ccl11","ccl13","ccl14","ccl15",
                                      "ccl16","ccl17","ccl18","ccl19",
                                      "ccl20","ccl21","ccl22","ccl23",
                                      "ccl24","ccl25","ccl26","ccl27",
                                      "ccl28","cxcl1","cxcl2","cxcl3",
                                      "cxcl4","cxcl4l1","cxcl5","cxcl6",
                                      "cxcl7","cxcl8","cxcl9","cxcl10",
                                      "cxcl11", "cxcl12", "cxcl13", "cxcl14", "cxcl16",
                                      "cxcl17", "cx3cl1", "xcl1", "xcl2")))
      levels(data$protein)
      data$protein <- factor(data$protein, levels = order.ck)
      
      data %>%
        #filter(pct_ortho >= 0.5) %>%
        ggplot() + 
        geom_tile(aes(protein, motif), fill = "steelblue3") +
        #scale_fill_gradient(low = "white", high = "steelblue3") +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
      # (5.2) CCL5 - 3-MER -----------------------------------------------------
      
      # import, reformat
      data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")  
      data <- data %>% filter(pct_ortho >= 0.5)
      
      # convert to matrix, fill NA with zero, convert back to long
      # data <- data %>% select(motif, protein, pct_ortho) %>% spread(protein, pct_ortho)
      # data[is.na(data)] <- 0
      # data <- data %>% gather(protein, pct_ortho, 2:ncol(data))
      
      # MASTER HEATMAP - MOTIF "SLICE"
      data <- data %>% filter(motif == "SPY" | motif == "PYS" | motif == "YSS" | 
                                motif == "SSD" | motif == "SDT" | motif == "DTT" | 
                                motif == "TTP" )
      
      order.mot <- as.factor(c("SPY","PYS","YSS","SSD","SDT","DTT","TTP"))
      levels(data$motif)
      data$motif <- factor(data$motif, levels = order.mot)
      
      # order chemokines
      order.ck <- as.factor(tolower(c("ccl5", "ccl1","ccl2", "ccl3", "ccl3l1",
                                      "ccl4","ccl4l1","ccl7","ccl8",
                                      "ccl11","ccl13","ccl14","ccl15",
                                      "ccl16","ccl17","ccl18","ccl19",
                                      "ccl20","ccl21","ccl22","ccl23",
                                      "ccl24","ccl25","ccl26","ccl27",
                                      "ccl28","cxcl1","cxcl2","cxcl3",
                                      "cxcl4","cxcl4l1","cxcl5","cxcl6",
                                      "cxcl7","cxcl8","cxcl9","cxcl10",
                                      "cxcl11", "cxcl12", "cxcl13", "cxcl14", "cxcl16",
                                      "cxcl17", "cx3cl1", "xcl1", "xcl2")))
      levels(data$protein)
      data$protein <- factor(data$protein, levels = order.ck)
      
      data %>%
        #filter(pct_ortho > 0.5) %>%
        ggplot() + 
        geom_tile(aes(protein, motif), fill = "steelblue3") +
        scale_fill_gradient(low = "white", high = ) +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
    # (5.3) ELR ----------------------------------------------------------------
      
      # import, reformat
      data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")  
      data <- subset(data, !grepl("C", data$motif)) # REMOVE CYSTEINE-CONTAINING MOTIFS
      raw <- read_csv("output/CK_UNSTRUCTURED_ALL_MERS.csv") # for reference
      
      # convert to matrix, fill NA with zero, convert back to long
      data <- data %>% select(motif, protein, pct_ortho) %>% spread(protein, pct_ortho)
      data[is.na(data)] <- 0
      data <- data %>% gather(protein, pct_ortho, 2:ncol(data))
      
      # MASTER HEATMAP - MOTIF "SLICE"
      data <- data %>% filter(motif == "EL" | motif == "LR" | motif == "ELR" )
      
      order.mot <- as.factor(c("EL","LR","ELR"))
      levels(data$motif)
      data$motif <- factor(data$motif, levels = order.mot)
      
      # order chemokines
      order.ck <- as.factor(tolower(c("ccl1","ccl2", "ccl3", "ccl3l1",
                                      "ccl4","ccl4l1","ccl5", "ccl7","ccl8",
                                      "ccl11","ccl13","ccl14","ccl15",
                                      "ccl16","ccl17","ccl18","ccl19",
                                      "ccl20","ccl21","ccl22","ccl23",
                                      "ccl24","ccl25","ccl26","ccl27",
                                      "ccl28","cxcl1","cxcl2","cxcl3",
                                      "cxcl4","cxcl4l1","cxcl5","cxcl6",
                                      "cxcl7","cxcl8","cxcl9","cxcl10",
                                      "cxcl11", "cxcl12", "cxcl13", "cxcl14", "cxcl16",
                                      "cxcl17", "cx3cl1", "xcl1", "xcl2")))
      levels(data$protein)
      data$protein <- factor(data$protein, levels = order.ck)
      
      data %>%
        #filter(pct_ortho > 0.5) %>%
        ggplot() + 
        geom_tile(aes(protein, motif, fill = pct_ortho)) +
        scale_fill_gradient(low = "white", high = "steelblue3") +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
     
      # (5.4) PxA / PxS --------------------------------------------------------
      
      # import, reformat
      data <- read_csv("output/CK_MOTIF_FREQUENCY.csv")  
      data <- data %>% filter(pct_ortho >= 0.5)
      
      # raw <- read_csv("output/CK_UNSTRUCTURED_ALL_MERS.csv") # for reference
      
      # convert to matrix, fill NA with zero, convert back to long
      data <- data %>% select(motif, protein, pct_ortho) %>% spread(protein, pct_ortho)
      data[is.na(data)] <- 0
      data <- data %>% gather(protein, pct_ortho, 2:ncol(data))
      
      # MASTER HEATMAP - MOTIF "SLICE"
      data <- data %>% filter(motif == "PxA" | motif == "PxS" | motif == "PxT" |
                                motif == "AA" | motif == "AS" | motif == "SA" | # "dummy motifs so that all chemokines are shown
                                motif == "AT" | motif == "EL" | motif == "LR" |
                                motif == "LS" | motif == "LL" | motif == "VA" |
                                motif == "SPY" | motif == "PYS" | motif == "YSS" | 
                                motif == "SSD" | motif == "SDT" | motif == "DTT" | 
                                motif == "TTP" |
                                motif == "SP" | motif == "PY" | motif == "YS" | 
                                motif == "SS" | motif == "SD" | motif == "DT" | 
                                motif == "TT" | motif == "TP" |
                                motif == "KPV" | motif == "PVS" | motif == "VSL" | 
                                motif == "SLS" | motif == "LSY" | motif == "SYR" |
                                motif == "AR" | motif == "AE" |motif == "SLS" |
                                motif == "ED" |motif == "SLS" | motif == "AK" |
                                motif == "GR" | motif == "FK" | motif == "LE" |
                                motif == "SK" | motif == "EG" | motif == "GV" )
      
      data <- data %>% mutate(yes.no = case_when(
        pct_ortho > 0 ~ 2,
        pct_ortho == 0 ~ 1
      ))
      
      order.mot <- as.factor(c("PxA","PxS", "PxT"))
      levels(data$motif)
      data$motif <- factor(data$motif, levels = rev(order.mot))
      
      # order chemokines
      order.ck <- as.factor(tolower(c("ccl1","ccl2", "ccl3", "ccl3l1",
                                      "ccl4","ccl4l1","ccl5", "ccl7","ccl8",
                                      "ccl11","ccl13","ccl14","ccl15",
                                      "ccl16","ccl17","ccl18","ccl19",
                                      "ccl20","ccl21","ccl22","ccl23",
                                      "ccl24","ccl25","ccl26","ccl27",
                                      "ccl28","cxcl1","cxcl2","cxcl3",
                                      "cxcl4","cxcl4l1","cxcl5","cxcl6",
                                      "cxcl7","cxcl8","cxcl9","cxcl10",
                                      "cxcl11", "cxcl12", "cxcl13", "cxcl14", "cxcl16",
                                      "cxcl17", "cx3cl1", "xcl1", "xcl2")))
      levels(data$protein)
      data$protein <- factor(data$protein, levels = order.ck)
      
      data %>%
        filter(pct_ortho > 0.5) %>%
        ggplot() + 
        geom_tile(aes(protein, motif), fill = "steelblue3") +
        #scale_fill_gradient(low = "white", high = "steelblue3") +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
      