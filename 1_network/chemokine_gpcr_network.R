# Name:     chemokine_gpcr_network.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 1F

##### PART 1: LOAD PACKAGES & SET WD ###########################################
  
  # load libraries; set workdir
  library(tidyverse)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/1_network/")


##### PART 2: FORMAT, WRITE TABLE ##############################################

  # matrix to long format
  data <- read_csv("input/network_matrix.csv")
  data <- data %>% 
    gather(ckr, score, 2:26) %>% 
    na.omit()
  write_csv(data, "output/network_matrix_long_raw.csv")    
  
  # Output table (output/network_matrix_long_raw.csv) was manually edited
  # to remove astrices and other characters for import into cytoscape
  # to make the final plot
  
  
  