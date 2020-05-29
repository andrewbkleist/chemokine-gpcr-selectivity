# Name:     5_conect_records.R
# Updated:  20191102
# User:     Andrew Kleist
# Figure:   Figure 6G

# This script genertes CONECT records that were manually imported into PyMol
# to make Figure 6G

##### LOAD PACKAGES & SET WD ###################################################
  
  # load libraries; set workdir
  library(tidyverse)
  require(Biostrings)
  library(ggrepel)
  library(bio3d)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_complex/")
  
##### FUNCTIONS ################################################################
  
  WriteCONECT <- function(RINFILE, PDBID, PDBFILE, OUTPUT){
    rin <- read_csv(RINFILE) %>% filter(file == PDBID) %>%
      mutate(type = case_when(
        Chain1 == "A" & Chain2 == "A" ~ "ck",
        Chain1 == "B" & Chain2 == "B" ~ "ckr",
        Chain1 == "A" & Chain2 == "B" ~ "inter"
      )) %>% filter(type =="inter")
    
    # read PDB, make df, select relevant columns
    pdb <- read.pdb(PDBFILE)
    pdb_df <- as.data.frame(pdb$atom)
    pdb_conv <- pdb_df %>% select(chain, resno, elety, eleno)
    pdb_conv <- pdb_conv %>% filter(elety == "CA")
    ck <- pdb_conv %>% filter(chain == "A")
    ckr <- pdb_conv %>% filter(chain == "B")
    
    
    # map atom indices to RIN file
    rin$ca1 <- ck$eleno[match(unlist(rin$ResNum1), ck$resno)]
    rin$ca2 <- ckr$eleno[match(unlist(rin$ResNum2), ckr$resno)]
    
    # clean up and write
    rin$CONECT <- c("CONECT")
    rin <- rin %>% select(CONECT, ca1, ca2)
    write_csv(rin, OUTPUT)
    
    # return
    return(rin)
    
    # remove
    rm(rin, pdb, pdb_df, ck, ckr, pdb_conv)
  }
  
##### 1: WRITE CONECT RECORDS ##################################################
  
  pdb.5uiw <- WriteCONECT("input/RIN.csv", 
                          "5uiw", 
                          "pdbs/5uiw_ck_clean.pdb", 
                          "conect/5uiw_conect.csv")
  
  pdb.4rws <- WriteCONECT("input/RIN.csv", 
                          "4rws", 
                          "pdbs/4rws_ck_clean.pdb", 
                          "conect/4rws_conect.csv")
  
  pdb.4xt1 <- WriteCONECT("input/RIN.csv", 
                          "4xt1", 
                          "pdbs/4xt1_ck_clean.pdb", 
                          "conect/4xt1_conect.csv")
  
  pdb.5wb2 <- WriteCONECT("input/RIN.csv", 
                          "5wb2", 
                          "pdbs/5wb2_clean.pdb", 
                          "conect/5wb2_conect.csv")
  
