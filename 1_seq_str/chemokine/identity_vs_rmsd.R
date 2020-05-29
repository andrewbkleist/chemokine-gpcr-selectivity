# Name:     identity_vs_rmsd.r
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 1B

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  library(RColorBrewer)
  library(bio3d)
  library(hexbin)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/1_seq_str/chemokine/")

##### ANALYSIS #################################################################
  
  # (1) PREPARE DATA FRAME -----------------------------------------------------
  
  # STRUCTURE
  # Import structure RMS info
  ck.struct <- read.table("input/ck_out_config.rms_rot")
  names <- ck.struct$V1
  ck.struct <- ck.struct[,-1]
  ck.struct$ck1 <- names
  colnames(ck.struct) <- names
  
  ck.struct <- gather(ck.struct, key = "ck2", value = "rms", 1:36)
  colnames(ck.struct) <- c("ck1", "ck2", "rms")
  
  # unite pairwise columns
  ck.struct <- unite_(ck.struct, "ck1ck2", c("ck1", "ck2"))
  
  
  # SEQUENCE
  # Make sequence x sequence identity
  # Ran unaligned, full (all but signal peptide) sequences through ClustalO:
  # https://www.ebi.ac.uk/Tools/msa/clustalo/
  # took output pairwise seq identities
  ck.seq <- read.table("input/identity_matrix_mod.txt", header = FALSE)
  names <- ck.seq[,1]
  rownames(ck.seq) <- ck.seq[,1]
  ck.seq <- ck.seq[,-1]
  colnames(ck.seq) <- rownames(ck.seq)
  
  # gather, unite
  ck.seq$ck1 <- names
  ck.seq <- gather(ck.seq, key = "ck2", value = "pid", 1:46)
  
  # unite pairwise columns
  ck.seq <- unite_(ck.seq, "ck1ck2", c("ck1", "ck2"))
  
  
  # COMBINE SEQUENCE AND STRUCTURE
  seq.struct <- inner_join(ck.struct, ck.seq)
  seq.struct$dist <- 100 - seq.struct$pid  
  seq.struct <- subset(seq.struct, seq.struct$pid != 100)
  seq.struct$ck1ck2 <- factor(seq.struct$ck1ck2, levels = seq.struct$ck1ck2[order(seq.struct$rms)])
  seq.struct$ck1ck2  # notice the changed order of factor levels
  

  # (2) PLOT DATA --------------------------------------------------------------
  
  # (1) 2D PLOT
  # Grey contour w dots
  # https://stackoverflow.com/questions/48282989/show-only-high-density-areas-with-ggplot2s-stat-density-2d
  ggplot(seq.struct, aes(rms, pid)) +
    geom_point(alpha=0.1) +
    ylim(0,100) +
    xlim(0,5) +
    stat_density_2d(geom = "polygon", aes(fill = ..level.., alpha = ..level..), bins = 30) +
    scale_fill_gradientn(colours = rev( brewer.pal( 7, "Greys" ) )) +
    scale_alpha_continuous(range = c(0.1, 0.3)) +
    theme_minimal()
  
  # (2) 1D PLOT - PID
  ggplot(seq.struct, aes(pid)) +
    geom_histogram(binwidth =2) +
    xlim(0,100) +
    coord_flip() +
    theme_minimal()
  
  # (2) 1D PLOT - RMSD
  ggplot(seq.struct, aes(rms)) +
    geom_histogram() +
    xlim(0,5) +
    coord_flip() +
    theme_minimal()
  
