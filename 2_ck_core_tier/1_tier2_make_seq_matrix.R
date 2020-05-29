# Name:     1_make_seq_matrix.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 2E 

# Note that this script makes the input for Figure 2E; does not make graph

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/2_ck_core_tier/")


##### 1: IMPORT SEQUENCES AND MAKE SEQUENCE MATRICES ###########################
  
  # 1.1: Import sequence, assign names
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
  
  ck <- Align2DataFrame("sequences/FULL_CHEMOKINE_ALIGNMENT_NODUPL.fasta",
                        "sequences/FULL_CHEMOKINE_CCN.txt")
  
  # 1.2: Select structured domains
  ck <- ck %>% select(seq, CX.1:H.10)
  
  # 1.3: Add chemokine, receptor labels
  pname <- ck %>% separate(seq, "/", remove = TRUE) %>% select(1)
  colnames(pname) <- c("protein")
  ck <- cbind(pname, ck)
  rm(pname)
  
  # 1.3: Add CC, CXC
  cc.cxc <- read_csv("input/cc_cxc_lookup.csv")
  ck$class <- cc.cxc$cc_cxc[match(unlist(ck$protein), cc.cxc$ckckr)]
  rm(cc.cxc)
  ck <- ck %>% select(protein, seq, class, CX.1:H.10)

  # 1.4: Filter CC, CXC
  ck <- ck %>% filter(class == "cc" | class == "cxc")

  # 1.5: Write CSV
  write_csv(ck, "output/chemokine_core_dataframe.csv")
