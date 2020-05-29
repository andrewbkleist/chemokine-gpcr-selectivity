# Name:     4_generate_chemokine_df.R
# Updated:  20191218
# Author:   Greg Slodkowicz / Andrew Kleist
# Figure:   Figure 6I,J (makes input)

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################
  
  library(tidyverse)
  library(ggrepel)
  library(ggforce)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")

##### 1: CHEMOKINE DATA ORGANIZATION ###########################################
# Creates DF with cancer, SNP, and Tier information (including motif vs. frag)

  # (1) IMPORT SNP, CANCER RAW COUNTS & COMBINE
  snp <- read_csv("output/CK_SNP_FREQ_MATRIX.csv")
  snp <- snp %>% gather(ccn, value, 2:ncol(snp))
  colnames(snp)[3] <- c("snp_count")
  
  cancer <- read_csv("output/CK_CANCER_MATRIX.csv")
  cancer <- cancer %>% gather(ccn, value, 2:ncol(cancer))
  colnames(cancer)[3] <- c("cancer_count")
  
  snp.cancer <- left_join(snp, cancer)
  rm(snp, cancer)
  
  # (2) ADD TIER INFORMATION
  # (2.1) Import and add tiers
  tiers <- read_csv("../2_ck_core_tier/output/MASTER_CONSERVATION_WITH_TIERS.csv") %>%
    dplyr::select(ck, ccn, tier, dom, ortho_cons, pctgaps)
  colnames(tiers)[1] <- c("protein")
  tiers$protein <- toupper(tiers$protein)
  
  snp.cancer <- left_join(snp.cancer, tiers)
  snp.cancer <- snp.cancer %>% gather(type, count, 3:4)
  
  # (2.2) No Tier information for  XCL1, XCL2, CX3CL1, CCL4L1
  snp.cancer <- snp.cancer %>% filter(protein != "XCL1" & 
                                        protein != "XCL2" &
                                        protein != "CX3CL1" &
                                        protein != "CCL4L1")  
  
  # (2.3) Add domain labels for NT and CT
  snp.cancer <- snp.cancer %>% mutate(dom = case_when(
    grepl('NT', ccn) ~ "NT",
    grepl('CT', ccn) ~ "CT",
    !grepl('NT', ccn) & !grepl('CT', ccn) ~ dom
  ))
  
  # (2.4) add "tier" label for NT
  snp.cancer <- snp.cancer %>% mutate(tier = case_when(
    dom == "NT" ~ "NT",
    dom != "NT" ~ tier
  ))  
  
  # (2.5) Add whether Tier 2 positions are conserved among orthologs
  snp.cancer <- snp.cancer %>% mutate(tier2_cons = case_when(
    tier == "tier2" & ortho_cons >= 0.8 ~ "yes",
    tier == "tier2" & ortho_cons < 0.8 ~ "no",
    tier != "tier2" ~ "non_tier2"
  ))
  
  # remove objects
  rm(tiers)
  
  # (3) REMOVE "NA" VALUES AND C-TERMINUS
  # NOTE THAT COUNTS CAN BE 0 OR NA -> "0" MEANS THAT THERE WAS A POSITION
  # THERE FOR WHICH NO SNPS/MUTATIONS OCCURED, WHEREAS "NA" MEANS THAT
  # THE CHEMOKINE DOES NOT HAVE A RESIDUE AT THAT POSITION
  snp.cancer <- snp.cancer %>% filter(!is.na(count))
  snp.cancer <- snp.cancer %>% filter(dom != "CT")
  
  # (4) ADD MOTIF INFORMATION
  # separate NT versus non
  nt <- snp.cancer %>% filter(dom == "NT")
  nonnt <- snp.cancer %>% filter(dom != "NT")
  
  # import motifs
  motif <- read_csv("output/ck_motif_conservation.csv") %>%
    dplyr::select(-motif, -pct_ortho) %>% unique()
  motif <- motif %>% mutate(value = case_when(
    tier ==  "motif" ~ 1,
    tier ==  "fragment" ~ 1
  ))
  motif <- motif %>% filter(!is.na(tier)) # removes XCL1, XCL2, CCL4L1
  motif <- motif %>% spread(tier, value)
  motif[is.na(motif)] <- 0
  motif <- motif %>% mutate(frag_or_no = case_when(
    fragment == 0 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 1 ~ "fragment",
    fragment == 0 & motif == 1 ~ "motif"
  ))
  motif <- motif %>% dplyr::select(-motif, -fragment)
  
  # join
  nt <- left_join(nt, motif)
  
  # remove positions where no "frag_or_no" ie, cleaved signal peptide positions
  nt <- nt %>% filter(!is.na(frag_or_no))
  nt <- nt %>% dplyr::select(-tier)
  
  nt$tier <- nt$frag_or_no
  nt <- nt %>% dplyr::select(-frag_or_no)
  
  colnames(nt)[4] <- c("ortho_cons")
  
  # recombine...
  snp.cancer <- bind_rows(nt, nonnt)
  rm(nt, nonnt, motif)
  
  
  # (5) ADD INTERFACE INFORMATION
  rin <- read_csv("../6_complex/input/RIN.csv") %>%
    mutate(type = case_when(
      Chain1 == "A" & Chain2 == "A" ~ "ck",
      Chain1 == "B" & Chain2 == "B" ~ "ckr",
      Chain1 == "A" & Chain2 == "B" ~ "inter"
    )) %>% filter(type =="inter") %>%
    filter(class == "xray") %>%
    dplyr::select(source_gnccn) %>% unique()
  rin <- rin$source_gnccn
  snp.cancer <- snp.cancer %>% mutate(interface = case_when(
    ccn %in% rin ~ "yes",
    !(ccn %in% rin)  ~ "no"
  ))
  
  # (6) ADD NUMBER OF ORTHOLOGS
  northo <- read_csv("../3_ck_motif/input/no_ortho.csv") %>% dplyr::select(-partners, -class)
  northo$ck <- toupper(northo$ck)
  colnames(northo)[1] <- c("protein")
  snp.cancer <- left_join(snp.cancer, northo)
  
  # (7) WRITE OUTPUT
  write_csv(snp.cancer, "output/CK_CANCER_SNP_TIER_OUTPUT.csv")
  rm(snp.cancer, rin, northo)
  