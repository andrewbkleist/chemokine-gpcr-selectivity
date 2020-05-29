# Name:     5_generate_receptor_df.R
# Updated:  20191218
# Author:   Greg Slodkowicz / Andrew Kleist
# Figure:   Figure 6 (makes input)

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################
  
  library(tidyverse)
  library(ggrepel)
  library(ggforce)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")

##### 1: RECEPTOR DATA ORGANIZATION ############################################
# Creates DF with cancer, SNP, and Tier information (including motif vs. frag)
  
  # (1) IMPORT SNP, CANCER RAW COUNTS & COMBINE
  snp <- read_csv("output/CKR_SNP_FREQ_MATRIX.csv")
  snp <- snp %>% gather(ccn, value, 2:ncol(snp))
  colnames(snp)[3] <- c("snp_count")
  
  cancer <- read_csv("output/CKR_CANCER_MATRIX.csv")
  cancer <- cancer %>% gather(ccn, value, 2:ncol(cancer))
  colnames(cancer)[3] <- c("cancer_count")
  
  snp.cancer <- left_join(snp, cancer)
  rm(snp, cancer)
  colnames(snp.cancer)[2] <- c("gn")
  
  # (2) ADD TIER INFORMATION
  # (2.1) Import and add tiers
  tiers <- read_csv("../4_ckr_core_tier/output/MASTER_CONSERVATION_WITH_TIERS.csv") %>%
    dplyr::select(ckr, gn, tier, dom, ortho_cons, pctgaps)
  colnames(tiers)[1] <- c("protein")
  tiers$protein <- toupper(tiers$protein)
  
  # make non-tier by combining low ortho and NA
  tiers <- tiers %>% mutate(tier = case_when(
    tier == "low_ortho" ~ "non_tier",
    is.na(tier) ~ "non_tier",
    tier != "low_ortho" & !is.na(tier) ~ tier
  ))
  
  # need to convert 1x22_ through 1x29_ nomenclature in "tiers" df to
  # NT.90 nomenclature...
  # recall that tiers are named in new nomenclature which encompasses
  # receptor N-terminal alignment including the conserved Cys, whereas
  # the SNP/cancer counts are counted with respect to the traditional GPCRdb
  # alignment; convert the TIER df to encompass the new nomenclature
  nt.conversion <- read_csv("input/lookup_for_new_gn_to_nt.csv") %>%
    filter(!is.na(old))
  
  tiers <- left_join(tiers, nt.conversion)
  tiers <- tiers %>% mutate(gn = case_when(
    is.na(old) ~ gn,
    !is.na(old) ~ old
  ))
  tiers <- tiers %>% dplyr::select(-old)
  rm(nt.conversion)
  # note that in HUMAN sequences, CXCR5 is the only one with 1x291-1x293
  # positions under the new nomenclature
  
  snp.cancer <- left_join(snp.cancer, tiers)
  snp.cancer <- snp.cancer %>% gather(type, count, 3:4)
  
  # (2.2) No Tier information for  CX3CR1, XCR1, ACKRs
  snp.cancer <- snp.cancer %>% filter(protein != "CX3CR1" & 
                                        protein != "CX3C1" &
                                        protein != "XXCR1" &
                                        protein != "XCR1" &
                                        protein != "ACKR1" &
                                        protein != "ACKR2" &
                                        protein != "ACKR3" &
                                        protein != "ACKR4" &
                                        protein != "CCRL2")
  
  
  # (2.3) Add domain labels for NT, CT, ECL1-3, ICL1-3
  # Also correct for NT nomenclature for positions now with TM1 designation
  snp.cancer <- snp.cancer %>% mutate(dom = case_when(
    grepl('NT', gn) & is.na(dom) & is.na(tier)  ~ "NT",
    grepl('CT', gn) ~ "CT",
    grepl('ECL1', gn) ~ "ECL1",
    grepl('ECL2', gn) ~ "ECL2",
    grepl('ECL3', gn) ~ "ECL3",
    grepl('ICL1', gn) ~ "ICL1",
    grepl('ICL2', gn) ~ "ICL2",
    grepl('ICL3', gn) ~ "ICL3",
    !grepl('NT', gn) & !grepl('CT', gn) & 
      !grepl('ECL1', gn) & !grepl('ECL2', gn) & !grepl('ECL3', gn) & 
      !grepl('ICL1', gn) & !grepl('ICL2', gn) & !grepl('ICL3', gn) ~ dom
  ))
  snp.cancer <- snp.cancer %>% mutate(dom = case_when(
    grepl('12x', gn) ~ "ICL1",
    grepl('23x', gn) ~ "ECL1",
    grepl('34x', gn) ~ "ICL2",
    grepl('45x', gn) ~ "ECL2",
    grepl('56x', gn) ~ "ICL3",
    grepl('67x', gn) ~ "ECL3",
    !grepl('12x', gn) & !grepl('23x', gn) & !grepl('34x', gn) & 
      !grepl('45x', gn) & !grepl('56x', gn) & !grepl('67x', gn) ~ dom
  ))
  snp.cancer <- snp.cancer %>% mutate(dom = case_when(
    grepl('gn1x', gn) ~ "TM1",
    !grepl('gn1x', gn) ~ dom
  ))
  
  # the only domains with NA are the NT regions that were renamed as TM1 segments
  # (ie including and after the conserved N-terminal Cys)
  snp.cancer <- snp.cancer %>% mutate(dom = case_when(
    is.na(dom) ~ "TM1",
    !is.na(dom) ~ dom
  ))
  
  # (2.6) Add whether Tier 2 positions are conserved among orthologs
  snp.cancer <- snp.cancer %>% mutate(tier2_cons = case_when(
    tier == "tier2" & ortho_cons >= 0.8 ~ "yes",
    tier == "tier2" & ortho_cons < 0.8 ~ "no",
    tier != "tier2" ~ "non_tier2"
  ))
  #snp.cancer <- snp.cancer %>% select(-ortho_cons)
  
  
  # remove objects
  rm(tiers)
  
  # (3) REMOVE "NA" VALUES AND C-TERMINUS
  # NOTE THAT COUNTS CAN BE 0 OR NA -> "0" MEANS THAT THERE WAS A POSITION
  # THERE FOR WHICH NO SNPS/MUTATIONS OCCURED, WHEREAS "NA" MEANS THAT
  # THE CHEMOKINE DOES NOT HAVE A RESIDUE AT THAT POSITION
  snp.cancer <- snp.cancer %>% filter(!is.na(count))
  # snp.cancer <- snp.cancer %>% filter(dom != "CT")
  
  # (4) ADD MOTIF INFORMATION
  
  # separate NT, ECL2 versus non
  nt <- snp.cancer %>% filter(dom == "NT")
  ecl2_unst <- snp.cancer %>% filter(grepl('gnECL2', gn))
  ecl2_str <- snp.cancer %>% filter(grepl('45x', gn))
  other <- snp.cancer %>% filter(dom != "NT" & !grepl('gnECL2', gn) & !grepl('45x', gn))
  other <- bind_rows(other, ecl2_str)
  rm(ecl2_str)  
  
  # import motifs NT
  motif.nt <- read_csv("output/ckr_nt_motif_conservation.csv") %>% 
    dplyr::select(-motif, -pct_ortho) %>% unique()
  motif.nt <- motif.nt %>% mutate(value = case_when(
    tier ==  "motif" ~ 1,
    tier ==  "fragment" ~ 1
  ))
  motif.nt <- motif.nt %>% spread(tier, value)
  motif.nt[is.na(motif.nt)] <- 0
  motif.nt <- motif.nt %>% mutate(frag_or_no = case_when(
    fragment == 0 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 1 ~ "fragment",
    fragment == 0 & motif == 1 ~ "motif"
  ))
  motif.nt <- motif.nt %>% select(-motif, -fragment)
  
  # import motifs ECL2
  motif.ecl2 <- read_csv("output/ckr_ecl2_motif_conservation.csv") %>% 
    dplyr::select(-motif, -pct_ortho) %>% unique()
  motif.ecl2 <- motif.ecl2 %>% mutate(value = case_when(
    tier ==  "motif" ~ 1,
    tier ==  "fragment" ~ 1
  ))
  motif.ecl2 <- motif.ecl2 %>% spread(tier, value)
  motif.ecl2[is.na(motif.ecl2)] <- 0
  motif.ecl2 <- motif.ecl2 %>% mutate(frag_or_no = case_when(
    fragment == 0 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 1 ~ "fragment",
    fragment == 0 & motif == 1 ~ "motif"
  ))
  motif.ecl2 <- motif.ecl2 %>% dplyr::select(-motif, -fragment)
  
  
  # add "gn" to front
  motif.nt$temp <- c("gn")
  motif.nt <- motif.nt %>% dplyr::select(protein, frag_or_no, temp, gn)
  motif.nt <- motif.nt %>% unite(new, 3:4, sep = "")
  colnames(motif.nt)[3] <- c("gn")
  
  motif.ecl2$temp <- c("gn")
  motif.ecl2 <- motif.ecl2 %>% dplyr::select(protein, frag_or_no, temp, gn)
  motif.ecl2 <- motif.ecl2 %>% unite(new, 3:4, sep = "")
  colnames(motif.ecl2)[3] <- c("gn")
  
  # join
  nt <- left_join(nt, motif.nt)
  ecl2 <- left_join(ecl2_unst, motif.ecl2) 
  # note that there are NA values in frag_or_no column of both nt and ecl2 dfs;
  # In NT these correspond to the motif preceding the cysteine (PC) which is not counted
  # Ditto for ECL2
  
  # remove positions where no "frag_or_no" position (see previous comment)
  nt <- nt %>% filter(!is.na(frag_or_no))
  nt <- nt %>% dplyr::select(-tier)
  nt$tier <- nt$frag_or_no
  nt <- nt %>% dplyr::select(-frag_or_no)
  
  ecl2 <- ecl2 %>% filter(!is.na(frag_or_no))
  ecl2 <- ecl2 %>% dplyr::select(-tier)
  ecl2$tier <- ecl2$frag_or_no
  ecl2 <- ecl2 %>% dplyr::select(-frag_or_no)
  
  colnames(nt)[4] <- c("pct_ortho")
  colnames(ecl2)[4] <- c("pct_ortho")
  
  # recombine...
  snp.cancer <- bind_rows(nt, ecl2, other)
  rm(nt, ecl2, motif.nt, motif.ecl2, ecl2_unst, other)
  
  # (4.5) ADD INTERFACE INFORMATION
  rin <- read_csv("../6_complex/input/RIN.csv") %>%
    mutate(type = case_when(
      Chain1 == "A" & Chain2 == "A" ~ "ck",
      Chain1 == "B" & Chain2 == "B" ~ "ckr",
      Chain1 == "A" & Chain2 == "B" ~ "inter"
    )) %>% filter(type =="inter") %>%
    filter(class == "xray") %>%
    dplyr::select(target_gnccn) %>% unique()
  #rin <- rin$target_gnccn
  
  # add gn
  rin$temp <- c("gn")
  rin <- rin %>% dplyr::select(temp, target_gnccn)
  rin <- rin %>% unite(new, 1:2, sep = "")
  colnames(rin) <- c("gn")
  
  # add 1x22, etc
  rin <- rin$gn
  rin <- c(rin, "gnNT.48", "gnNT.49", "gnNT.50", "gnNT.51", "gnNT.52")
  
  snp.cancer <- snp.cancer %>% mutate(interface = case_when(
    (gn %in% rin) & dom == "TM1" ~ "yes",
    (gn %in% rin) & dom == "NT" ~ "no",
    (gn %in% rin) & dom != "TM1" ~ "yes",
    !(gn %in% rin)  ~ "no"
  ))
  
  
  # ADD non_tier label
  snp.cancer <- snp.cancer %>% mutate(tier = case_when(
    !is.na(tier) ~ tier,
    is.na(tier) ~ "non_tier"
  ))
  
  # (6) WRITE OUTPUT
  write_csv(snp.cancer, "output/CKR_CANCER_SNP_TIER_OUTPUT.csv")
  rm(snp.cancer, rin)
