# Name:     get_gpcr_ligands.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 1E


##### LOAD PACKAGES & SET WD ###################################################

  # load libraries; set workdir
  library(tidyverse)
  library(ggrepel)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/1_iuphar/")


##### ANALYSIS #################################################################
  
  # import data
  data <- read_csv("input/gpcr_ligand.csv") %>% 
    select(-rat_gene_symbol, -mouse_gene_symbol, -family)
  
  # count numbers of ligands and receptors per family
  ligands <- data %>% select(family2, ligand, class) %>% distinct() %>% count(family2)  
  receptors <- data %>% select(iuphar, family2, class) %>% distinct() %>% count(family2) 
  colnames(ligands)[2] <- c("n_ligands")
  colnames(receptors)[2] <- c("n_receptors")

  # join numbers of ligands and receptors to a single df
  final <- left_join(ligands, receptors) 
  final$class <- data$class[match(unlist(final$family2), data$family2)]
  final$ratio <- log2(final$n_ligands / final$n_receptors)
  
  # plot
  final %>%
    ggplot(aes(n_receptors, n_ligands, color=class)) +
    geom_point() +
    xlim(0,50) +
    ylim(0,50) +
    geom_point(shape = 21, color = "grey20", fill = "white",size = 2, stroke = 0.5) +
    geom_text_repel(data=subset(final, n_receptors>5 & n_ligands>5), aes(label=family2)) +
    theme_minimal()
    theme_minimal()
    