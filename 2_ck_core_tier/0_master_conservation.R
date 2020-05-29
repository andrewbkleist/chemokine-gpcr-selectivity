# Name:     0_master_conservation.R
# Updated:  20191031
# User:     Andrew Kleist
# Figure:   Figure 2D

# Note that this script generates conservation scores for each position in the 
# master alignment of all 1056 sequences and identifies Tier 1 positions

##### LOAD PACKAGES & SET WD ###################################################
 
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ape)
  library(stringr)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/2_ck_core_tier/")

  
##### PART 1: CONSERVATION SCORING FULL ALIGNMENT ##############################
  
  # (1.1) RUN CONSERVATION -----------------------------------------------------
  
  # load fasta files, convert to matrix, assign CCN
  chemokine_master <- readAAMultipleAlignment("sequences/FULL_CHEMOKINE_ALIGNMENT.fasta")
  lookup_master <- as.matrix(chemokine_master)

  ccn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
  colnames(lookup_master) <- c(ccn_names)

  # define alignment input and output for conservation scoring, run scoring
  ck_names <- c("FULL_CHEMOKINE_ALIGNMENT.fasta")
  ck_names2 <- c("ALL_ortho_para.txt")
  
  setwd("sequences/")
  # call mstatx from command line (MstatX must be installed locally and 
  # the appropriate path to the scoring matrix must be specified)
  j <- 1 # create loop index for ck_names3 assignment (ck_names will use 'i')
  for(i in ck_names){
    tryit <- paste('mstatx -i ', print(i), '-s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ', print(ck_names2[j]))
    system(tryit)
    j <- j + 1
  }
  setwd("..")
  
  
  # (1.2) IMPORT AND TIDY ------------------------------------------------------
  
  # import trident scores to R  
  setwd("sequences/")
  fetchfile <- Sys.glob("ALL_ortho_para.txt")
  trident <- lapply(fetchfile, read.table)
  names(trident) <- Sys.glob("ALL_ortho_para.txt")
  rm(fetchfile)
  setwd("..")
  
  # create matrix of all trident scores (row) by chemokine (column)
  tritemp <- sapply(trident, "[[",2)
  trimat <- as.data.frame(tritemp)
  ccn_names <- as.data.frame(t(t(ccn_names)))
  colnames(ccn_names) <- c("ccn")
  rownames(ccn_names) <- 1:nrow(ccn_names)
  trimat <- cbind(ccn_names, trimat)
  rm(tritemp, trident, ck_names, ck_names2, i, j, tryit)
  ck.scores <- trimat
  colnames(ck.scores)[2] <- c("value")
  
  # (1.3) ADD ADDITIONAL DESCRIPTORS -------------------------------------------
  
  # import attribute files
  attr.ck <- read_csv("input/bwccn_domain_conversion.csv")
  colnames(attr.ck)[1] <- c("ccn")
  
  # change ccn to character
  ck.scores$ccn <- as.character(ck.scores$ccn)
  
  # merge with existing
  ck.scores <- left_join(ck.scores, attr.ck)
  
  # count gap numbers, calculate percent gaps per position
  gaps <- NULL
  for (i in 1:ncol(lookup_master)){
    temp <- as.data.frame(table(lookup_master[,i]))
    colnames(temp) <- c("resno", "n")
    temp <- subset(temp, resno == "-")
    temp$pos <- i
    gaps <- rbind.data.frame(gaps, temp)
    rm(temp)
  }
  rm(i)
  gaps$pct <- gaps$n/nrow(lookup_master)
  gaps <- cbind(ccn_names, gaps)
  colnames(gaps) <- c("ccn", "restype", "n", "pos", "pct")
  gaps <- select(gaps, ccn, n, pct)
  gaps$ccn <- as.character(gaps$ccn)
  
  # merge with existing
  ck.scores <- left_join(ck.scores, gaps)
  
  # reorder, colnames
  colnames(ck.scores)[4] <- c("ngaps")
  colnames(ck.scores)[5] <- c("pctgaps")
  
  
  # write output file
  write_csv(ck.scores, "output/MASTER_CONSERVATION_ALL_1058.csv")
  
  rm(attr.ck, ccn_names, chemokine_master, ck.scores, gaps, lookup_master,
       name, trident, trimat)

      
##### PART 2: IDENTIFY TIER 1 ##################################################
  
  data <- read_csv("output/MASTER_CONSERVATION_ALL_1058.csv")
  data <- data[96:189,] # select CCN positions in core only (CX.1-H.10)

  # add Tier 1 designation
  data <- data %>% 
    mutate(tier = case_when(
      data$value >= 0.8  ~ "tier_1")) 

  # order and plot
  data$ccn <- factor(data$ccn, levels = data$ccn[order(data$value)])
  data %>% filter(value >= 0.5) %>%
    ggplot(aes(ccn, value)) +
    geom_point() +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()  
    
  