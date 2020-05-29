# Name:     4_assign_tiers.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 4C


##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ape)
  library(stringr)
  library(ggrepel)
  library(seqinr)
  library(bio3d)
  library(seqRFLP)
  library(ggpubr)
  library(RColorBrewer)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/2_ck_core_tier/")

  
##### PART 1: IDENTIFY AND PLOT TIER 1 & TIER 2 ################################
  
  # import data, etc
  data <- read.csv("output/CHEMOKINE_CONSERVATION.csv")
  data <- data %>% filter(dom != "NTc" & dom != "CT")   # remove NT, CT
  
  # add class information
  class <- read_csv("input/cc_cxc_lookup.csv")
  data$class <- class$cc_cxc[match(unlist(data$ck), class$ckckr)]
  rm(class)
  
  # (1.1): ADD CLASSIFICATION LABELS ---------------------------------------------
  
  # remove non CC, CXC
  data <- data %>% filter(class == "cc" | class == "cxc")
  
  # Get Tier 1
  tier1 <- c("CX.1", "CX.5", "b1b2.12", "B3.3", "b3h.2", "H.3")
  
  # Get Tier 2
  tier2 <- read_csv("output/CHEMOKINE_CLASSIFICATION_N3.csv") %>%
    filter(mean >= 0.8) # %>% filter(int == "int")
  tier2 <- tier2$motif
  
  # add labels
  data <- data %>%
    mutate(tier = case_when(
      data$ccn %in% tier1 ~ "tier1",
      data$ccn %in% tier2 ~ "tier2",
      data$cc < 0.8 & data$ortho_cons >= 0.8 & data$class == "cc" ~ "tier3",
      data$cxc < 0.8 & data$ortho_cons >= 0.8 & data$class == "cxc" ~ "tier3",
      data$ortho_cons < 0.8  ~ "low_ortho"
    ))
  
  write_csv(data, "output/MASTER_CONSERVATION_WITH_TIERS.csv")
  
  # (1.2): ORDER AND PLOT --------------------------------------------------------
  
  # order by CCN, chemokine
  order.ck <- as.factor(rev(tolower(c("CCL1","CCL2", "CCL3", "CCL3L1",
                                  "CCL4","CCL4L1","CCL5","CCL7","CCL8",
                                  "CCL11","CCL13","CCL14","CCL15",
                                  "CCL16","CCL17","CCL18","CCL19",
                                  "CCL20","CCL21","CCL22","CCL23",
                                  "CCL24","CCL25","CCL26","CCL27",
                                  "CCL28","CXCL1","CXCL2","CXCL3",
                                  "CXCL4","CXCL4L1","CXCL5","CXCL6",
                                  "CXCL7","CXCL8","CXCL9","CXCL10",
                                  "CXCL11","CXCL12","CXCL13", "CXCL14", "CXCL16",
                                  "CXCL17", "CX3CL1", "XCL1", "XCL2"))))
  levels(data$ck)
  data$ck <- factor(data$ck, levels = order.ck)
  
  order.ccn <- as.factor(c("CX.1","CX.2","CX.3","CX.4","CX.5","cxb1.1","cxb1.2","cxb1.3","cxb1.4","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.9","cxb1.10","cxb1.11","cxb1.12","cxb1.13","cxb1.14","cxb1.15","cxb1.16","cxb1.17","cxb1.18","cxb1.19","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.1","b1b2.2","b1b2.3","b1b2.4","b1b2.5","b1b2.6","b1b2.7","b1b2.8","b1b2.9","b1b2.10","b1b2.11","b1b2.12","b1b2.13","b1b2.14","b1b2.15","b1b2.16","b1b2.17","b1b2.18","b1b2.19","b1b2.20","b1b2.21","b1b2.22","b1b2.23","b1b2.24","b1b2.25","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.1","b2b3.2","b2b3.3","b2b3.4","b2b3.5","b2b3.6","b2b3.7","b2b3.8","b2b3.9","b2b3.10","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","b3h.5","b3h.6","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10"))
  levels(data$ccn)
  data$ccn <- factor(data$ccn, levels = order.ccn)
  
  
  # plot CORE
  data %>%
    filter(pctgaps < 0.6) %>% # only plotting positions for which gaps <60%
    filter(ck != "cx3cl1") %>%
    ggplot() + 
    geom_tile(aes(ccn, ck, fill = tier))+
    scale_fill_manual(values=c("grey90", "gray40", "steelblue4", "mediumslateblue")) +
    coord_fixed() +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  # write positions for alignment plot
  positions <- data %>% filter(ck == "ccl11") %>% filter(pctgaps < 0.6) 
    # arbitrary choice of CCL11
  positions <- positions$ccn

  
  # (1.3): COUNT NUMBER T3 FOR PLOT ----------------------------------------------
  
  # count no. Tier 3
  no.tier <- data %>% dplyr::select(ck, tier) #%>% count(ck, tier)
  no.tier <- no.tier %>% dplyr::count(ck, tier)
  
  # count no Tier 3 at interface
  rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
    filter(Chain1 != Chain2) %>% select(source_gnccn) %>% distinct()
  inter <- rin$source_gnccn

  data <- data %>% mutate(interface = case_when(
    data$ccn %in% inter ~ "yes"
  ))
  data$interface[is.na(data$interface)] <- c("no")
  
  # now count
  no.tier <- data %>% dplyr::select(ck, tier, interface) #%>% count(ck, tier)
  no.tier <- no.tier %>% dplyr::count(ck, tier, interface)
  
  # (1.4): SUBSET MASTER ALIGNMENT AND WRITE FASTA FOR FIGURE ------------------
  chemokine_master <- readAAMultipleAlignment("sequences/FULL_CHEMOKINE_ALIGNMENT.fasta")
  lookup_master <- as.matrix(chemokine_master)
  lookup_master <- as.data.frame(lookup_master)
  
  ccn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
  colnames(lookup_master) <- c(ccn_names)
  
  aln <- lookup_master[(names(lookup_master) %in% positions)]
  aln <- aln[1:43,]
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
  
  writeFasta(aln.unite, "output/example.fasta")
  

  