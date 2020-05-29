# Name:     3_make_conservation.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 2; Supp. Figure 4C

# Note that this script makes the input for Figure and associated Supp. figures; 
# does not make graph

##### LOAD PACKAGES & SET WD ###################################################
 
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ape)
  library(stringr)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/2_ck_core_tier/")

  
##### PART 1: GENERATE ALIGNMENTS FOR SCORING ##################################
# Generate orthologous sequence sets from the master alignment
  
  # load fasta files, convert to matrix, assign CCN
  chemokine_master <- readAAMultipleAlignment("sequences/FULL_CHEMOKINE_ALIGNMENT.fasta")
  lookup_master <- as.matrix(chemokine_master)

  ccn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
  colnames(lookup_master) <- c(ccn_names)

  ck_list <- c("human",
               "ccl1 ", "ccl2 ", "ccl3 ", "ccl3l1 ", "ccl4 ", "ccl5 ",
               "ccl7 ", "ccl8 ", "ccl11 ", "ccl13 ", "ccl14 ", "ccl15 ", "ccl16 ", 
               "ccl17 ", "ccl18 ", "ccl19 ", "ccl20 ", "ccl21 ", "ccl22 ", "ccl23 ",
               "ccl24 ", "ccl25 ", "ccl26 ", "ccl27 ", "ccl28 ", "cxcl1 ", 
               "cxcl2 ", "cxcl3 ", "cxcl4 ", "cxcl4l1 ", "cxcl5 ", "cxcl6 ", 
               "cxcl7 ", "cxcl8 ", "cxcl9 ", "cxcl10 ", "cxcl11 ", "cxcl12 ", 
               "cxcl13 ", "cxcl14 ", "cxcl16 ", "cxcl17 ", "cx3cl1 ")
    
    # 43 chemokines (missing XCL1, XCL2), 44 entries total including "ALL_para"
  
  ck_names <- c("ALL_para.fasta",
                "ccl1_ortho.fasta", "ccl2_ortho.fasta", "ccl3_ortho.fasta",
                "ccl3l1_ortho.fasta","ccl4_ortho.fasta", "ccl5_ortho.fasta",
                "ccl7_ortho.fasta", "ccl8_ortho.fasta", "ccl11_ortho.fasta",
                "ccl13_ortho.fasta", "ccl14_ortho.fasta", "ccl15_ortho.fasta",
                "ccl16_ortho.fasta", "ccl17_ortho.fasta", "ccl18_ortho.fasta",
                "ccl19_ortho.fasta", "ccl20_ortho.fasta", "ccl21_ortho.fasta",
                "ccl22_ortho.fasta", "ccl23_ortho.fasta", "ccl24_ortho.fasta",
                "ccl25_ortho.fasta", "ccl26_ortho.fasta", "ccl27_ortho.fasta",
                "ccl28_ortho.fasta", "cxcl1_ortho.fasta", "cxcl2_ortho.fasta",
                "cxcl3_ortho.fasta", "cxcl4_ortho.fasta", "cxcl4l1_ortho.fasta",
                "cxcl5_ortho.fasta", "cxcl6_ortho.fasta", "cxcl7_ortho.fasta",
                "cxcl8_ortho.fasta", "cxcl9_ortho.fasta", "cxcl10_ortho.fasta",
                "cxcl11_ortho.fasta", "cxcl12_ortho.fasta", "cxcl13_ortho.fasta",
                "cxcl14_ortho.fasta", "cxcl16_ortho.fasta", "cxcl17_ortho.fasta",
                "cx3cl1_ortho.fasta")
  
  # loop to write fasta files
  j <- 1 # create loop index for ck_names assignment (ck_list will use 'i')
  setwd("sequences/")
  for(i in ck_list){
    ck_rows <- grep(i, rownames(chemokine_master))
    ck_head <- head(ck_rows, n = 1)
    ck_tail <- tail(ck_rows, n = 1)
    rowmask(chemokine_master, invert = TRUE) <- IRanges(start = ck_head, end = ck_tail) 
    # mask non-selected rows
    temp_aln <- as(chemokine_master, "AAStringSet") # cast as AA string set
    writeXStringSet(temp_aln, file = ck_names[j]) # write fasta
    rowmask(chemokine_master) <- NULL # unmask chemokine_master alignment
    j <- j + 1
  }

  # remove objects
  rm(i, j, ck_head, ck_list, ck_names, ck_rows, ck_tail, temp_aln, ccn_names, chemokine_master, lookup_master)
  
  setwd("..")
  
##### PART 2: CONSERVATION SCORING #############################################
# From all alignments, calculate per-position conservation scores
# Note that CC and CXC sequence sets were manually curated and not generated
# by the loop in the section above
  
  # define sequence alignments to score
  ck_names <- c("ALL_para.fasta", "CC_para.fasta", "CXC_para.fasta", 
                "ccl1_ortho.fasta", "ccl2_ortho.fasta", "ccl3_ortho.fasta",
                "ccl3l1_ortho.fasta","ccl4_ortho.fasta", "ccl5_ortho.fasta",
                "ccl7_ortho.fasta", "ccl8_ortho.fasta", "ccl11_ortho.fasta",
                "ccl13_ortho.fasta", "ccl14_ortho.fasta", "ccl15_ortho.fasta",
                "ccl16_ortho.fasta", "ccl17_ortho.fasta", "ccl18_ortho.fasta",
                "ccl19_ortho.fasta", "ccl20_ortho.fasta", "ccl21_ortho.fasta",
                "ccl22_ortho.fasta", "ccl23_ortho.fasta", "ccl24_ortho.fasta",
                "ccl25_ortho.fasta", "ccl26_ortho.fasta", "ccl27_ortho.fasta",
                "ccl28_ortho.fasta", "cxcl1_ortho.fasta", "cxcl2_ortho.fasta",
                "cxcl3_ortho.fasta", "cxcl4_ortho.fasta", "cxcl4l1_ortho.fasta",
                "cxcl5_ortho.fasta", "cxcl6_ortho.fasta", "cxcl7_ortho.fasta",
                "cxcl8_ortho.fasta", "cxcl9_ortho.fasta", "cxcl10_ortho.fasta",
                "cxcl11_ortho.fasta", "cxcl12_ortho.fasta", "cxcl13_ortho.fasta",
                "cxcl14_ortho.fasta", "cxcl16_ortho.fasta", "cxcl17_ortho.fasta",
                "cx3cl1_ortho.fasta")  
  
  # define sequence alignments to score
  ck_names2 <- c("ALL_para.txt", "CC_para.txt", "CXC_para.txt", 
                 "ccl1_ortho.txt", "ccl2_ortho.txt", "ccl3_ortho.txt",
                 "ccl3l1_ortho.txt","ccl4_ortho.txt", "ccl5_ortho.txt",
                 "ccl7_ortho.txt", "ccl8_ortho.txt", "ccl11_ortho.txt",
                 "ccl13_ortho.txt", "ccl14_ortho.txt", "ccl15_ortho.txt",
                 "ccl16_ortho.txt", "ccl17_ortho.txt", "ccl18_ortho.txt",
                 "ccl19_ortho.txt", "ccl20_ortho.txt", "ccl21_ortho.txt",
                 "ccl22_ortho.txt", "ccl23_ortho.txt", "ccl24_ortho.txt",
                 "ccl25_ortho.txt", "ccl26_ortho.txt", "ccl27_ortho.txt",
                 "ccl28_ortho.txt", "cxcl1_ortho.txt", "cxcl2_ortho.txt",
                 "cxcl3_ortho.txt", "cxcl4_ortho.txt", "cxcl4l1_ortho.txt",
                 "cxcl5_ortho.txt", "cxcl6_ortho.txt", "cxcl7_ortho.txt",
                 "cxcl8_ortho.txt", "cxcl9_ortho.txt", "cxcl10_ortho.txt",
                 "cxcl11_ortho.txt", "cxcl12_ortho.txt", "cxcl13_ortho.txt",
                 "cxcl14_ortho.txt", "cxcl16_ortho.txt", "cxcl17_ortho.txt",
                 "cx3cl1_ortho.txt")
  
  setwd("sequences/")
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
  
##### PART 3: IMPORT AND TIDY SCORES ###########################################
# Create master data frame with all conservation scores
  
  # (3.1) IMPORT AND TIDY ------------------------------------------------------
  
    # import trident scores to R
    setwd("sequences/")
    fetchfile <- Sys.glob("*.txt")
    trident <- lapply(fetchfile, read.table)
    names(trident) <- Sys.glob("*.txt")
    rm(fetchfile)
    setwd("..")
    
    # import sequence, ccn
    chemokine_master <- readAAMultipleAlignment("sequences/FULL_CHEMOKINE_ALIGNMENT.fasta")
    lookup_master <- as.matrix(chemokine_master)
    
    ccn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
    colnames(lookup_master) <- c(ccn_names)
    
    # create matrix of all trident scores (row) by chemokine (column)
    tritemp <- sapply(trident, "[[",2)
    trimat <- as.data.frame(tritemp)
    trimat <- trimat[, -5] # remove CCN column (artifact from global *txt import)
    ccn_names <- as.data.frame(t(t(ccn_names)))
    colnames(ccn_names) <- c("ccn")
    rownames(ccn_names) <- 1:nrow(ccn_names)
    trimat <- cbind(ccn_names, trimat)
    rm(tritemp)
    
    # tidy data
    ck.scores <- gather(trimat, key = "ck", value = "ortho_cons", 6:ncol(trimat))
    colnames(ck.scores) <- c("ccn", "all_ortho", "all_para", "cc", "cxc", "file", "ortho_cons")
    ck.scores <- select(ck.scores, file, ccn, all_ortho, all_para, cc, cxc, ortho_cons) 
    
  # (3.2) ADD ADDITIONAL DESCRIPTORS -------------------------------------------
  
    # import attribute files
    attr.ck <- read_csv("input/bwccn_domain_conversion.csv")
    colnames(attr.ck)[1] <- c("ccn")
    
    # change ccn to character
    ck.scores$ccn <- as.character(ck.scores$ccn)
    
    # merge with existing
    ck.scores <- left_join(ck.scores, attr.ck)

    # add number of gaps per position and calculate gap pct
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
    ck.scores <- ck.scores %>% select(file,ccn,dom,all_ortho,all_para,cc,cxc,ortho_cons,n,pct)
    colnames(ck.scores)[9] <- c("ngaps")
    colnames(ck.scores)[10] <- c("pctgaps")
    
    # change to chemokine name from file name
    name <- strsplit(ck.scores$file, "_", fixed = TRUE)
    name <- t(as.data.frame(lapply(name, "[", 1)))
    rownames(name) <- 1:nrow(name)
    colnames(name) <- c("ck")
    
    # bind
    ck.scores <- cbind(name, ck.scores)
    
    # write output file
    write_csv(ck.scores, "output/CHEMOKINE_CONSERVATION.csv")
  
    rm(attr.ck, ccn_names, chemokine_master, ck.scores, gaps, lookup_master,
       name, trident, trimat)
    
