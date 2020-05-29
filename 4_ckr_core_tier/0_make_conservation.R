# Name:     0_make_conservation.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 4; Supp. Figure 7; Supp. Figure 8

# Note that this script performs conservation scoring to make data used
# in Figures 4, Supp. Figures 7 and 8, but does not make graph

##### LOAD PACKAGES & SET WD ###################################################
 
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ape)
  library(stringr)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")

  
##### PART 1: CONSERVATION SCORING #############################################
# From all alignments, calculate per-position conservation scores
# Note that all sequence sets were manually curated
  
  # (1.1) RUN CONSERVATION -----------------------------------------------------
  
  # define file names
  ck_names <- c("ALL_para.fasta", "CC_para.fasta", "CXC_para.fasta", "ACK_para.fasta", 
                "FULL_RECEPTOR_ALIGNMENT_CYS_ADJ.fasta",
                "ccr1_ortho.fasta", "ccr2_ortho.fasta", "ccr3_ortho.fasta",
                "ccr4_ortho.fasta","ccr5_ortho.fasta", "ccr6_ortho.fasta",
                "ccr7_ortho.fasta", "ccr8_ortho.fasta", "ccr9_ortho.fasta",
                "ccr10_ortho.fasta", "cxcr1_ortho.fasta", "cxcr2_ortho.fasta",
                "cxcr3_ortho.fasta", "cxcr4_ortho.fasta", "cxcr5_ortho.fasta",
                "cxcr6_ortho.fasta", "cx3cr1_ortho.fasta", "xcr1_ortho.fasta",
                "ackr1_ortho.fasta", "ackr2_ortho.fasta", "ackr3_ortho.fasta",
                "ackr4_ortho.fasta", "ccrl2_ortho.fasta")  
  
  ck_names2 <- c("ALL_para.txt", "CC_para.txt", "CXC_para.txt", "ACK_para.txt", 
                 "FULL_RECEPTOR_ALIGNMENT_CYS_ADJ.txt",
                 "ccr1_ortho.txt", "ccr2_ortho.txt", "ccr3_ortho.txt",
                 "ccr4_ortho.txt","ccr5_ortho.txt", "ccr6_ortho.txt",
                 "ccr7_ortho.txt", "ccr8_ortho.txt", "ccr9_ortho.txt",
                 "ccr10_ortho.txt", "cxcr1_ortho.txt", "cxcr2_ortho.txt",
                 "cxcr3_ortho.txt", "cxcr4_ortho.txt", "cxcr5_ortho.txt",
                 "cxcr6_ortho.txt", "cx3cr1_ortho.txt", "xcr1_ortho.txt",
                 "ackr1_ortho.txt", "ackr2_ortho.txt", "ackr3_ortho.txt",
                 "ackr4_ortho.txt", "ccrl2_ortho.txt")  
  
  setwd("sequences_for_mstatx/")
  # call mstatx from command line
  j <- 1 # create loop index for ck_names3 assignment (ck_names will use 'i')
  for(i in ck_names){
    tryit <- paste('mstatx -i ', print(i), '-s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ', print(ck_names2[j]))
    system(tryit)
    j <- j + 1
  }
  setwd("..")
  
  # remove variables
  rm(ck_names, ck_names2, i, j, tryit)
  
  
##### PART 2: IMPORT AND TIDY SCORES ###########################################
  
  # (2.1) IMPORT AND TIDY ------------------------------------------------------
  
  # import trident scores to R
  setwd("sequences_for_mstatx/")
  fetchfile <- Sys.glob("*.txt")
  trident <- lapply(fetchfile, read.table)
  names(trident) <- Sys.glob("*.txt")
  rm(fetchfile)
  setwd("..")
  
  # import sequence, ccn
  chemokine_master <- readAAMultipleAlignment("sequences_for_mstatx/FULL_RECEPTOR_ALIGNMENT_CYS_ADJ.fasta")
  lookup_master <- as.matrix(chemokine_master)
  
  ccn_names <- c(read.table("sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt", sep = ",", colClasses = "character"))
  colnames(lookup_master) <- c(ccn_names)
  
  # create matrix of all trident scores (row) by chemokine (column)
  tritemp <- sapply(trident, "[[",2)
  trimat <- as.data.frame(tritemp)
  ccn_names <- as.data.frame(t(t(ccn_names)))
  colnames(ccn_names) <- c("gn")
  rownames(ccn_names) <- 1:nrow(ccn_names)
  trimat <- cbind(ccn_names, trimat)
  rm(tritemp)
  
  # tidy data
  ck.scores <- gather(trimat, key = "ckr", value = "ortho_cons", 7:ncol(trimat))
  colnames(ck.scores) <- c("gn", "ack", "all_para", "cc", "cxc", "all_951", "file", "ortho_cons")
  ck.scores <- select(ck.scores, file, gn, all_951, all_para, cc, cxc, ack, ortho_cons) 
  
  # (1.3) ADD ADDITIONAL DESCRIPTORS -------------------------------------------
  
  # import attribute files
  attr.ck <- read_csv("input/bwccn_domain_conversion.csv")
  colnames(attr.ck)[1] <- c("gn")
  
  # change ccn to character
  ck.scores$gn <- as.character(ck.scores$gn)
  
  # merge with existing
  ck.scores <- left_join(ck.scores, attr.ck)
  
  # add gap numbers and percentages
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
  colnames(gaps) <- c("gn", "restype", "n", "pos", "pct")
  gaps <- select(gaps, gn, n, pct)
  gaps$gn <- as.character(gaps$gn)
  
  # merge with existing
  ck.scores <- left_join(ck.scores, gaps)
  
  # reorder, colnames
  ck.scores <- ck.scores %>% select(file,gn,dom,all_para, all_951,cc,cxc,ack,ortho_cons,n,pct)
  colnames(ck.scores)[10] <- c("ngaps")
  colnames(ck.scores)[11] <- c("pctgaps")
  
  # change to chemokine name from file name
  name <- strsplit(ck.scores$file, "_", fixed = TRUE)
  name <- t(as.data.frame(lapply(name, "[", 1)))
  rownames(name) <- 1:nrow(name)
  colnames(name) <- c("ckr")
  
  # bind
  ck.scores <- cbind(name, ck.scores)
  
  # write output file
  write_csv(ck.scores, "output/RECEPTOR_CONSERVATION.csv")
  
  rm(attr.ck, ccn_names, ck.scores, gaps, lookup_master,
     name, trident, trimat, chemokine_master)
  
  