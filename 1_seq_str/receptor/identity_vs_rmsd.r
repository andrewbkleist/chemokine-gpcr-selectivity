# Name:     identity_vs_rmsd.r
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 1C

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  library(RColorBrewer)
  library(bio3d)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/1_seq_str/receptor/")

##### ANALYSIS #################################################################
  
  # (1) PREPARE DATA FRAME -----------------------------------------------------
  
  # STRUCTURE
  # Import structure RMS info
  ckr.struct <- read.csv("pdbs/ckr_out_config.rms_rot_20191220.csv")
  ckr.struct <- gather(ckr.struct, key = "ckr2", value = "rms", 2:6)
  ckr.struct <- unite_(ckr.struct, "ckr1ckr2", c("ckr1", "ckr2"))
  
  # SEQUENCE
  # Make sequence x sequence identity
  # Ran unaligned core sequences through ClustalO:
  # https://www.ebi.ac.uk/Tools/msa/clustalo/
  # took output pairwise seq identities
  ckr.seq <- read.csv("input/clustal_pid_matrix.csv", header = FALSE)
  names <- ckr.seq$V1
  rownames(ckr.seq) <- ckr.seq[,1]
  ckr.seq <- ckr.seq[,-1]
  colnames(ckr.seq) <- rownames(ckr.seq)

  # gather, unite
  ckr.seq$ckr1 <- names
  ckr.seq <- gather(ckr.seq, key = "ckr2", value = "pid", 1:4)
  ckr.seq <- unite_(ckr.seq, "ckr1ckr2", c("ckr1", "ckr2"))
  
  # COMBINE SEQUENCE AND STRUCTURE
  seq.struct <- inner_join(ckr.struct, ckr.seq)
  seq.struct$dist <- 100 - seq.struct$pid  

  # (2) PLOT DATA --------------------------------------------------------------
  
  # PLOT - 2D
  ggplot(seq.struct, aes(rms, pid)) +
    geom_point(aes(alpha = 0.8)) +
    ylim(0,100) +
    xlim(0,5) +
    #stat_ellipse() +
    #geom_hex(bins = 10) +
    geom_text(aes(label=ckr1ckr2),hjust=0, vjust=0) +
    theme_minimal()

  # PLOT - 1D RMSD
  ggplot(seq.struct, aes(rms)) +
    geom_histogram() +
    xlim(0,5) +
    #coord_flip() +
    #geom_vline(xintercept=68.46043) +
    theme_minimal()
  
##### PART 2: ALL CHEMOKINE RECEPTOR SEQUENCE HISTOGRAM ########################
  
  # SEQUENCE
  # Make sequence x sequence identity
  # Ran unaligned core sequences through ClustalO:
  # https://www.ebi.ac.uk/Tools/msa/clustalo/
  # took output pairwise seq identities
  all.ckr.seq <- read.table("input/ckr_clustal_pid_matrix.csv", header = FALSE)
  names <- all.ckr.seq[,1]
  rownames(all.ckr.seq) <- all.ckr.seq[,1]
  all.ckr.seq <- all.ckr.seq[,-1]
  colnames(all.ckr.seq) <- rownames(all.ckr.seq)
  
  # gather, unite
  all.ckr.seq$ckr1 <- names
  all.ckr.seq <- gather(all.ckr.seq, key = "ckr2", value = "pid", 1:23)
  all.ckr.seq <- unite_(all.ckr.seq, "ckr1ckr2", c("ckr1", "ckr2"))
  
  all.ckr.seq$dist <- 100 - all.ckr.seq$pid  
  all.ckr.seq <- subset(all.ckr.seq, all.ckr.seq$pid != 100)
  
  # PLOT
  ggplot(all.ckr.seq, aes(pid)) +
    geom_histogram(binwidth =2) +
    xlim(0,100) +
    coord_flip() +
    #geom_vline(xintercept=68.46043) +
    theme_minimal()
 