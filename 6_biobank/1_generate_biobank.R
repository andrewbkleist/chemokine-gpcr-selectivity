# Name:     1_generate_snp.R
# Updated:  20200317
# Author:   Greg Slodkowicz, Andrew Kleist
# Figure:     Figure 6 (generates input only)

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################

  library(tidyverse)
  library(plyr)
  library(seqinr)
  library(stringr)
  library(gplots)
  library(Biostrings)
  AMINO_ACID_CODE <- toupper(AMINO_ACID_CODE)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_biobank/")
  
##### FUNCTIONS ################################################################
  
  # FUNCTION 1 -----------------------------------------------------------------
  format_variants_gnomad <- function(var_table) {
    
    # define variants
    variant_types <- c("missense", "frameshift", "inframe insertion", "inframe deletion", "stop gained", "stop lost", "start lost")
    variant_colours <- c("blue", "indianred1",  "indianred2", "yellow", "red", "cyan2", "cyan3")
    names(variant_colours) <- variant_types
    
    var_subset <- 
      subset(var_table, !Annotation %in% c("5' UTR", "3' UTR", "synonymous", "intron",
                                           "non coding transcript exon", 
                                           "splice donor", "splice region",
                                           "splice acceptor", "upstream gene", 
                                           "mature miRNA", "downstream gene",
                                           "stop retained"))
    
    # filter low confidence
    var_subset$Flags[is.na(var_subset$Flags)] <- ""
    var_subset <- subset(var_subset, substr(Flags, 1, 2) != "LC")
    var_subset <- subset(var_subset, !(Annotation == "frameshift" & Protein.Consequence == ""))
    var_subset$MAF <- var_subset$Allele.Frequency
    var_subset$MAF[var_subset$MAF > 0.5] <- 1-var_subset$MAF[var_subset$MAF > 0.5]
    
    var_subset$MAF_cat <- factor(">=5%", levels=c("<1/10,000", "<1/1,000", "<1%", "<5%", ">=5%"))
    var_subset$MAF_cat[var_subset$MAF < 5/100] <- "<5%"
    var_subset$MAF_cat[var_subset$MAF < 1/100] <- "<1%"
    var_subset$MAF_cat[var_subset$MAF < 1/1000] <- "<1/1,000"
    var_subset$MAF_cat[var_subset$MAF < 1/10000] <- "<1/10,000"
    print("MAF categories assigned")
    
    # make the consequence names more pithy
    var_subset$Annotation[var_subset$Annotation == "missense_variant"] <- "missense"
    var_subset$Annotation[var_subset$Annotation == "frameshift_variant"] <- "frameshift"
    
    var_subset$Annotation <- factor(var_subset$Annotation, levels=variant_types)
    aa_consequence <- 
      str_match(var_subset$Protein.Consequence, "p.([A-Z][a-z][a-z])([0-9]+)(.+)")[, 2:4]
    aa_consequence <- 
      data.frame(aa_ref=aa_consequence[, 1], aa_pos=aa_consequence[, 2], aa_alt=aa_consequence[, 3], stringsAsFactors=F)
    var_subset[, c("aa_ref", "aa_pos", "aa_alt")] <- aa_consequence
    var_subset$aa_pos <- as.numeric(var_subset$aa_pos)
    
    var_subset
  }
  
  # FUNCTION 2 -----------------------------------------------------------------
  get_aligned_variants <- function(r) {
    gene_symbol <- r$gene_symbol
    print(gene_symbol)
    gnomad_file <- Sys.glob(paste0("gnomad/", gene_symbol, "_*.csv"))[1]
    print(gnomad_file)
    gnomad_vcf <- read.csv(gnomad_file, stringsAsFactors=F, na.strings="NA")
    gnomad_subset <- format_variants_gnomad(gnomad_vcf)
    # print(tail(gnomad_subset))
    
    aln_split <- str_split(r$aln, "", simplify=T)
    gaps <- aln_split == "-"
    seq <- paste0(aln_split[!gaps])
    gap_offsets <- rep(NA, length(aln_split))
    gap_offsets[!gaps] <- 1:length(seq)
    
    gnomad_subset <- subset(gnomad_subset, Annotation == "missense")
    print(subset(gnomad_subset, aa_pos > length(seq)))
    cbind(gnomad_subset, 
          data.frame(gene_symbol=gene_symbol, 
                     aln_pos=sapply(gnomad_subset$aa_pos, function(i) { which(i == gap_offsets) }), stringsAsFactors=F))
  }
  
  # FUNCTION 3 -----------------------------------------------------------------
  make_variant_matrix <- function(aln, variants, allele.count=F, threeletter=T) {
    
    # (1) LOOP #1 - MAKE EMPTY MATRIX FOR MASTER CHEMOKINE ALIGNMENT
    # acts on "aln" which corresponds to ccl_seqs
    variant_matrix <- matrix(0, nrow=length(aln), ncol=nchar(aln[1])) # make blank matrix
    rownames(variant_matrix) <- names(aln) # change rownames
    
    for (s in rownames(variant_matrix)) {  # for each row...
      seq <- aln[s]                        # define new var that is a vector of 
      # seq positions from that row
      
      for (i in 1:(nchar(aln[1]))) {       # for each seq position...
        if (substr(seq, i, i) == '-') {    # replace dashes (ie no AA in alignment)
          variant_matrix[s, i] <- NA       # with NA
        }
      }
    }                                      # after this loop, all matrix positions
    # with zeros indicate that a chemokine
    # has a residue there, and all NAs
    # indicate that a chemokine has a dash
    
    # (2) LOOP #2 - MULTIPLE MANIPULATIONS OF EMPTY (ZERO) MATRIX
    # acts on "aln" which corresponds to ccl_seqs AND acts on ccl_variants which
    # is the master SNP table from GNOMAD
    
    for(i in 1:nrow(variants)) {          # for each row in GNOMAD table
      r <- variants[i, ]                  # make vector of values from table
      seq <- aln[r$gene_symbol]           # grab gene name (eg CCL1) and seq...
      
      
      # (2.1) IF-ELSE STATEMENT - IF "threeletter" IS TRUE
      # In both cases, fetches amino acid - outputs "K" for instance
      if (threeletter)
        aa_aln <- AMINO_ACID_CODE[toupper(substr(seq, r$aln_pos, r$aln_pos))] 
      else
        aa_aln <- toupper(substr(seq, r$aln_pos, r$aln_pos))
      
      
      # (2.2) IF STATEMENT - FIND MISMATCHES IN LISTED SNP AND REFERENCE
      if(toupper(r$aa_ref) != aa_aln) {
        print(paste("Mismatch in", names(seq), "at", r$aa_pos, "Gnomad:", toupper(r$aa_ref), "Aln:", aa_aln))
      }
      
      
      # (2.3) IF-ELSE STATEMENT - IF "allele.count" IS TRUE
      if (allele.count) # if "allele.count" is TRUE...
        variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + r$Allele.Count
      # replace the matrix of 0's at a specific chemokine (row) AND 
      # alignment position (column) WITH existing value (ie 0) plus 
      # MAF (defined in function 1)
      else
        variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + 1
    }
    # otherwise add 1 (indicating frequency of 100% unless a SNP was listed)
    variant_matrix
  }
  
  # FUNCTION 4 -----------------------------------------------------------------
  get_aligned_biobank <- function(r) {
    gene_symbol <- r$gene_symbol
    print(gene_symbol)
    biobank_subset <- subset(rbind(ccl_biobank, ccr_biobank), Gene == r$gene_symbol)
    if(nrow(biobank_subset) == 0)
      return (data.frame())
    
    aln_split <- str_split(r$aln, "", simplify=T)
    gaps <- aln_split == "-"
    seq <- paste0(aln_split[!gaps])
    gap_offsets <- rep(NA, length(aln_split))
    gap_offsets[!gaps] <- 1:length(seq)
    
    print(subset(biobank_subset, aa_pos > length(seq)))
    
    cbind(biobank_subset, 
          data.frame(gene_symbol=gene_symbol, 
                     aln_pos=sapply(biobank_subset$aa_pos, function(i) { which(i == gap_offsets) }), stringsAsFactors=F))
  }
  
  
##### 1: SNP TABLE AND MATRIX GENERATION #######################################
  
  # (1.1) IMPORT SEQS, MAKE SEQUENCE VECTORS ----------------------------------- 
  
  # CHEMOKINE - read sequences - in this case Nterm, Cterm, core independent fasta files 
  aln.nterm <- read.fasta("sequences/chemokine_paralogs_NTERM.fasta", as.string=T)
  aln.core <- read.fasta("sequences/chemokine_paralogs_CORE.fasta", as.string=T)
  aln.cterm <- read.fasta("sequences/chemokine_paralogs_CTERM.fasta", as.string=T)
  aln.nterm <- aln.nterm[names(aln.core)]
  aln.cterm <- aln.cterm[names(aln.core)]
  
  # create seq character strings AND write fasta
  ccl_seqs <- sapply(names(aln.core), function(n) { paste0(aln.nterm[[n]], aln.core[[n]], aln.cterm[[n]]) })
  write.fasta(as.list(ccl_seqs), names=names(aln.core), file.out="sequences/chemokine_paralogs.fasta")
  
  # RECEPTOR - read sequences
  aln.ccr <- read.fasta("sequences/chemokine_receptor_paralogs.fasta", as.string=T)
  ccr_seqs <- sapply(aln.ccr, as.character)
  
  # remove used objects
  rm(aln.cterm, aln.nterm)
  
  
  # (1.2) CHEMOKINE - PARSE GNOMAD - RAW DATA ----------------------------------
  
  # get chemokine GENE NAMES
  gene_symbols <- toupper(str_split(names(aln.core), "_", simplify=T)[, 1]) # get chemokine names
  names(ccl_seqs) <- gene_symbols # set chemokine names
  ccl_table <- data.frame(gene_symbol=gene_symbols, aln=ccl_seqs, stringsAsFactors=F) # create gene name / sequence table
  
  # removed used objects
  rm(gene_symbols, aln.core) # chemokine objects
  
  
  # (1.3) RECEPTOR - PARSE GNOMAD - RAW DATA ----------------------------------
  
  # get receptor GENE NAMES
  ccr_gene_symbols <- toupper(str_split(names(aln.ccr), "_", simplify=T)[, 1])
  names(ccr_seqs) <- ccr_gene_symbols
  ccr_table <- data.frame(gene_symbol=ccr_gene_symbols, aln=ccr_seqs, stringsAsFactors=F)
  
  # removed used objects
  rm(aln.ccr, ccr_gene_symbols) # receptor objects
  rm(AMINO_ACID_CODE, format_variants_gnomad, get_aligned_variants)
  
  
  # (1.4) BIOBANK MAPPINGS -----------------------------------------------------

  # Get a mapping of gene symbols to transcripts from Gnomad file names (elegant, I know)
  id_map <- str_split(dir(path="gnomad", pattern = "*.csv"), "_", simplify=T)
  id_map <- data.frame(gene_symbol=id_map[, 1], transcript_id=id_map[, 2], stringsAsFactors=F) %>%
    filter(gene_symbol %in% c(rownames(ccl_table), rownames(ccr_table)))
  
  ccl_biobank <- read.csv("biobank/chemokines_GA_hits_1e-2_processed_split.csv", sep=";", stringsAsFactors=F, header=T)
  ccr_biobank <- read.csv("biobank/receptors_GA_hits_1e-2_processed_split.csv", sep=";", stringsAsFactors=F, header=T)
  
  ccl_biobank <- ccl_biobank %>%
    mutate(Gene=toupper(Gene)) %>%
    filter( (Transcript_ID %in% id_map$transcript_id)) %>%
    mutate(aa_ref=substr(AA_Consequences, 1, 1),
           aa_alt=substr(AA_Consequences, nchar(AA_Consequences), nchar(AA_Consequences)),
           aa_pos=as.numeric(substr(AA_Consequences, 2, nchar(AA_Consequences)-1)))
  
  ccr_biobank <- ccr_biobank %>%
    mutate(Gene=toupper(Gene)) %>%
    filter(Transcript_ID %in% id_map$transcript_id) %>%
    mutate(aa_ref=substr(AA_Consequences, 1, 1),
           aa_alt=substr(AA_Consequences, nchar(AA_Consequences), nchar(AA_Consequences)),
           aa_pos=as.numeric(substr(AA_Consequences, 2, nchar(AA_Consequences)-1)))
  
  ccl_biobank_aln <- adply(ccl_table, 1, get_aligned_biobank)
  ccr_biobank_aln <- adply(ccr_table, 1, get_aligned_biobank)
  
  
  # (1.5) BIOBANK CCN/GN, TIDY -------------------------------------------------
  
  # filter for significiant
  ccl_biobank_aln <- ccl_biobank_aln %>% filter(P_value <= 1e-8)
  ccr_biobank_aln <- ccr_biobank_aln %>% filter(P_value <= 1e-8)
  
  # import GN CCN positions
  ccn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
  ccn_names <- as.data.frame(ccn_names)
  ccn_names <- t(ccn_names)
  colnames(ccn_names) <- c("gnccn")
  ccn_names <- as.data.frame(ccn_names)
  ccn_names$aln_pos <- 1:nrow(ccn_names)
  
  gn_names <- c(read.table("sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt", sep = ",", colClasses = "character"))
  gn_names <- as.data.frame(gn_names)
  gn_names <- t(gn_names)
  colnames(gn_names) <- c("gnccn")
  gn_names <- as.data.frame(gn_names)
  gn_names$aln_pos <- 1:nrow(gn_names)
  
  # add new column to biobank data
  ccl_biobank_aln$gnccn <- ccn_names$gnccn[match(unlist(ccl_biobank_aln$aln_pos), ccn_names$aln_pos)]
  ccr_biobank_aln$gnccn <- gn_names$gnccn[match(unlist(ccr_biobank_aln$aln_pos), gn_names$aln_pos)]
  rm(gn_names, ccn_names)
  rm(ccl_biobank, ccr_biobank)
  
  # remove CT positions
  ccl_biobank_aln <- ccl_biobank_aln %>% filter(!grepl('CT', gnccn) )
  ccr_biobank_aln <- ccr_biobank_aln %>% filter(!grepl('CT', gnccn) )
  
  # remove chemokines/receptors for which no Tier 2, 3 positions exist
  # (none for chemokines)
  ccr_biobank_aln <- ccr_biobank_aln %>% 
    filter(!grepl('ACKR', gene_symbol) ) %>%
    filter(gene_symbol != "CCRL2") %>%
    filter(gene_symbol != "XCR1")
    
  
  # (1.6) ADD TIER INFORMATION -------------------------------------------------
  
  # import tier, nontier, filter
  ck.tier <- read_csv("output/CK_CANCER_SNP_TIER_OUTPUT.csv")
  ck.tier <- ck.tier %>% select(protein, ccn, dom, tier, interface) %>% unique()
  colnames(ck.tier)[1] <- c("Gene")
  colnames(ck.tier)[2] <- c("gnccn")
  
  ckr.tier <- read_csv("output/CKR_CANCER_SNP_TIER_OUTPUT.csv")
  ckr.tier <- ckr.tier %>% select(protein, gn, dom, tier, interface) %>% unique()
  colnames(ckr.tier)[1] <- c("Gene")
  colnames(ckr.tier)[2] <- c("gnccn")
  
  # add tier, nontier
  ccl_biobank_aln <- left_join(ccl_biobank_aln, ck.tier)
  ccr_biobank_aln <- left_join(ccr_biobank_aln, ckr.tier)
  rm(ck.tier, ckr.tier)
  
  # remove nontier and fragment
  ccl_biobank_aln <- ccl_biobank_aln %>% 
    filter(tier != "low_ortho") %>%
    filter(tier != "fragment") %>%
    filter(!is.na(tier) )
  
  ccr_biobank_aln <- ccr_biobank_aln %>% 
    filter(tier != "non_tier") %>%
    filter(tier != "fragment") %>%
    filter(!is.na(tier) )
  
  # # make matrix
  # ccl_biobank_matrix <- make_variant_matrix(ccl_seqs, ccl_biobank_aln, allele.count=F, threeletter=F)
  # ccr_biobank_matrix <-make_variant_matrix(ccr_seqs, ccr_biobank_aln, allele.count=F, threeletter=F)
  
  rm(ccl_table, ccr_table, id_map, ccl_seqs, ccr_seqs, get_aligned_biobank, make_variant_matrix)

  # tidy and combine for table
  ccl_biobank_aln$type <- c("chemokine")
  ccr_biobank_aln$type <- c("receptor")
  
  master <- bind_rows(ccl_biobank_aln, ccr_biobank_aln)
  master <- master %>% select(-aln, -gene_symbol)
  write_csv(master, "output/Supplementary_Table_2.csv")
    