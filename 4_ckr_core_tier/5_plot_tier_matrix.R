# Name:     5_plot_tier_matrix.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 8A


##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ggrepel)
  library(bio3d)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")

##### 1: PLOT TIER MATRIX ######################################################
  
  # import data
  data <- read_csv("output/MASTER_CONSERVATION_WITH_TIERS.csv")
  
  # import interface
  rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
    filter(Chain1 != Chain2) %>% select(target_gnccn) %>% distinct()
  rin$gn <- c("gn")
  inter <- unite(rin, gn, sep = "", c(gn, target_gnccn))
  inter <- inter$gn
  inter <- c(inter,  "gn1x28_", "gn1x24_", "gn1x25_")
  rm(rin)
  
  # order receptors
  order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                   "ccr5","ccr6","ccr7","ccr8","ccr9",
                                   "ccr10","cxcr1","cxcr2","cxcr3",
                                   "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                   "ackr1","ackr2","ackr3","ackr4",
                                   "ccrl2")))
  data$ckr <- factor(data$ckr, levels = rev(order.ckr))
  
  # order positions
  order <- as.factor(unique(data$gn))
  data$gn <- factor(data$gn, levels = (order))
  
  # PLOT INTERFACE
  data %>%
    #filter(pctgaps < 0.6) %>%
    filter(gn %in% inter) %>%
    filter(ckr != "cx3cr1") %>%
    filter(ckr != "xcr1") %>%
    #filter(gn != "gn45x50" & gn != "gn45x51" & gn != "gn45x52") %>%
    #na.omit(tier) %>%
    ggplot() + 
    geom_tile(aes(gn, ckr, fill = tier))+
    scale_fill_manual(values=c("grey90", "gray40", "steelblue4", "mediumslateblue")) +
    # scale_fill_manual(values=c( "gray40", "steelblue4", "darkcyan", "dodgerblue1","aquamarine2","white")) +
    coord_fixed() +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
##### 2: WRITE CORRESPONDING SEQUENCES #########################################
  
  inter2 <- data %>% filter(gn %in% inter) %>% select(gn) %>% distinct()
  inter2 <- as.character(inter2$gn)
  
  chemokine_master <- readAAMultipleAlignment("sequences/FULL_RECEPTOR_ALIGNMENT_CYS_ADJ.fasta")
  lookup_master <- as.matrix(chemokine_master)
  lookup_master <- as.data.frame(lookup_master)
  
  ccn_names <- c(read.table("sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt", sep = ",", colClasses = "character"))
  colnames(lookup_master) <- c(ccn_names)
  
  aln <- lookup_master[(names(lookup_master) %in% inter2)] 
  # note all 21 non gn interface positions are dropped (eg gnNTr.Cp1)
  aln <- aln[grepl("*human", rownames(aln)), ]
  
  aln.unite <- unite(aln[,1:ncol(aln)], col = seq,  sep = "")
  aln.unite$name <- rownames(aln.unite)
  aln.unite <- select(aln.unite, name, seq)
  rownames(aln.unite) <- 1:nrow(aln.unite)
  
  writeFasta<-function(data, filename){
    fastaLines = c()
    for (rowNum in 1:nrow(data)){
      fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
      fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
    }
    fileConn<-file(filename)
    writeLines(fastaLines, fileConn)
    close(fileConn)
  }
  
  writeFasta(aln.unite, "output/RECEPTOR_INTERFACE_PARALOGS.fasta")
  