# Name:     general_ck_features.r
# Updated:  20191031
# User:     Andrew Kleist
# Figure:   Figure 9B

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  library(Peptides)
  library(seqinr)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/5_isoelectric/")

##### 1: IMPORT DATA ###########################################################
  
  # load fasta files, run analysis
  ck.seq <- read.fasta("ALL_para_core.fasta",  seqtype = c("AA"))
      # Get seq: ck.seq$`sp|P22362|24-96|CCL1`
  ck.stat <- lapply(ck.seq, AAstat)
      # ck.stat[[1]][1]  accesses number of each res
      # ck.stat[[1]][2]  accesses summary stats
  
  # compile stats into df: No. RESID
  ck.res.summ <- NULL
  for(i in 1:length(ck.stat)){
    temp <- as.data.frame(ck.stat[[i]][1])
    temp$ck <- names(ck.stat[i])
    colnames(temp) <- c("resid", "no", "ck")
    ck.res.summ <- rbind(ck.res.summ, temp)
    rm(temp)
  }
  rm(i)

  # compile stats into df: SUMMARY STATS
  ck.prop.summ <- NULL
  for(i in 1:length(ck.stat)){
    temp <- as.data.frame(ck.stat[[i]][2])
    temp <- gather(temp, stat, value, 1:9)
    temp$ck <- names(ck.stat[i])
    ck.prop.summ <- rbind(ck.prop.summ, temp)
    rm(temp)
  }
  rm(i)
  
  # compile stats into df: PI
  ck.pi.summ <- NULL
  for(i in 1:length(ck.stat)){
    temp <- as.data.frame(ck.stat[[i]][3])
    temp$ck <- names(ck.stat[i])
    ck.pi.summ <- rbind(ck.pi.summ, temp)
    rm(temp)
  }
  rm(i)
  
  ck.mw <- as.data.frame(mw(ck.seq))
  ck.mw$ck <- rownames(ck.mw)
  colnames(ck.mw) <- c("mw", "ck")
  
 
  
##### 1: PLOT DATA #############################################################
  
  # single stat plotting PI
  ck.pi.summ$type <- c("ck")
  ggplot(ck.pi.summ, aes(x = type, y = Pi)) +
    #geom_bar(stat = "identity") +
    geom_violin(trim = FALSE) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
    #coord_flip() +
    theme_minimal()
    
  