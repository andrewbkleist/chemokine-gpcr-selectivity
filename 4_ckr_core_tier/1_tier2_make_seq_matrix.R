# Name:     1_tier2_make_seq_matrix.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 4B; Supp. Figure 7

# Note that this script performs conservation scoring to make data used
# in Figures 4B, Supp. Figures 7 but does not make graph

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")


##### 1: IMPORT SEQUENCES AND MAKE SEQUENCE MATRICES ###########################
  
  # (1.1): Import sequence, assign names
  Align2DataFrame <- function(ALIGNMENT, NAMES){
    aln <- readAAMultipleAlignment(ALIGNMENT)
    aln.df <- as.matrix(aln)
    seqname <- as.tibble(rownames(aln.df))
    colnames(seqname) <- c("seq")
    aln.df <- as.tibble(aln.df)
    gnccn <- t(read.table(NAMES, sep = ","))
    colnames(aln.df) <- gnccn
    aln.df <- cbind(seqname, aln.df)
    
    return(aln.df)
    rm(aln, aln.df, gnccn, seqname)
  }
  
  ckr <- Align2DataFrame("sequences/FULL_RECEPTOR_ALIGNMENT_CYS_ADJ.fasta", 
                         "sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt")
  
  
  # (1.2): Select structured domains
  ckr <- ckr %>% select(seq,
                        gnNTr.Cm10:gn1x60, gn12x48:gn12x51, # TM1, ICL1
                        gn2x36:gn2x66, gn23x49:gn23x52, # TM2, ECL1
                        gn3x21:gn3x56, gn34x50:gn34x58, # TM3, ICL2
                        gn4x35:gn4x65, gn45x50:gn45x52, # TM4, ECL2
                        gn5x31:gn5x68, # TM5
                        gn6x28:gn6x65, # TM6
                        gn7x23:gn7x56) # TM7
  
  pname <- ckr %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  ckr <- cbind(pname, ckr)
  rm(pname)
  
  # (1.3): Add CC, CXC
  cc.cxc <- read_csv("input/cc_cxc_lookup.csv")
  ckr$class <- cc.cxc$cc_cxc[match(unlist(ckr$protein), cc.cxc$ckckr)]
  rm(cc.cxc)
  ckr <- ckr %>% select(protein, seq, class, gnNTr.Cm10:gn7x56)
  
  # (1.4): Filter CC, CXC
  ckr.cc.cxc <- ckr %>% filter(class == "cc" | class == "cxc")

  # (1.5): Write CSV
  write_csv(ckr.cc.cxc, "output/receptor_core_dataframe.csv")
  write_csv(ckr, "output/receptor_core_dataframe_ack.csv")
  
  # Name:     1_tier2_make_seq_matrix.R
  # Updated:  20191030
  # User:     Andrew Kleist
  # Figure:   Figure 4B; Supp. Figure 7
  
  # Note that this script performs conservation scoring to make data used
  # in Figures 4B, Supp. Figures 7 but does not make graph
  
 
##### 2: IMPORT ALIGNMENT WITH B2AR, ETC, FOR TIER 2 TESTING AGAINST NON-CKR ###
  
  # (2.1): Import sequence, assign names
  Align2DataFrame <- function(ALIGNMENT, NAMES){
    aln <- readAAMultipleAlignment(ALIGNMENT)
    aln.df <- as.matrix(aln)
    seqname <- as.tibble(rownames(aln.df))
    colnames(seqname) <- c("seq")
    aln.df <- as.tibble(aln.df)
    gnccn <- t(read.table(NAMES, sep = ","))
    colnames(aln.df) <- gnccn
    aln.df <- cbind(seqname, aln.df)
    
    return(aln.df)
    rm(aln, aln.df, gnccn, seqname)
  }
  
  ckr <- Align2DataFrame("sequences/FULL_RECEPTOR_ALIGNMENT_CYS_ADJ_B2AR_ETC.fasta", 
                         "sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt")
  
  
  # (2.2): Select structured domains
  ckr <- ckr %>% select(seq,
                        gnNTr.Cm10:gn1x60, gn12x48:gn12x51, # TM1, ICL1
                        gn2x36:gn2x66, gn23x49:gn23x52, # TM2, ECL1
                        gn3x21:gn3x56, gn34x50:gn34x58, # TM3, ICL2
                        gn4x35:gn4x65, gn45x50:gn45x52, # TM4, ECL2
                        gn5x31:gn5x68, # TM5
                        gn6x28:gn6x65, # TM6
                        gn7x23:gn7x56) # TM7
  
  pname <- ckr %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  ckr <- cbind(pname, ckr)
  rm(pname)
  
  # (2.3): Add CC, CXC
  cc.cxc <- read_csv("input/cc_cxc_lookup.csv")
  ckr$class <- cc.cxc$cc_cxc[match(unlist(ckr$protein), cc.cxc$ckckr)]
  rm(cc.cxc)
  ckr <- ckr %>% select(protein, seq, class, gnNTr.Cm10:gn7x56)
  
  ckr <- ckr %>% mutate(class = case_when(
    is.na(class) ~ "non",
    !is.na(class) ~ class
  ))

  write_csv(ckr, "output/receptor_core_dataframe_ack_b2ar.csv")
  
  