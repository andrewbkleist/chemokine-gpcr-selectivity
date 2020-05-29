# Name:     2_generate_cancer.R
# Updated:  20191218
# Author:   Greg Slodkowicz, Andrew Kleist
# Figure:     Figure 6 (generates input only)

##### LOAD PACKAGES, LIBRARIES, SET WD #########################################
  
  library(tidyverse)
  library(plyr)
  library(seqinr)
  library(stringr)
  library(gplots)
  library(Biostrings)
  library(biomaRt)
  AMINO_ACID_CODE <- toupper(AMINO_ACID_CODE)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")
  

##### FUNCTIONS ################################################################

  # FUNCTION 1 -----------------------------------------------------------------
  idx_codon <- function(idx) {
  # go from index to codon number and position within
    
    codon_n <- floor((idx-1) / 3) + 1
    codon_pos <- ((idx-1) %% 3) + 1
    
    return(c(codon_n, codon_pos))
  }
  
  # FUNCTION 2 -----------------------------------------------------------------
  get_aligned_cancer_variants <- function(r) {
    gene_symbol <- r$gene_symbol
    print(gene_symbol)
    cancer_muts_subset <- subset(cancer_muts_parsed, gene_symbol == r$gene_symbol)
    
    if (!nrow(cancer_muts_subset)) {
      return(data.frame())
    }
    
    aln_split <- str_split(r$aln, "", simplify=T)
    gaps <- aln_split == "-"
    seq <- paste0(aln_split[!gaps])
    gap_offsets <- rep(NA, length(aln_split))
    gap_offsets[!gaps] <- 1:length(seq)
    
    print(subset(cancer_muts_subset, aa_pos > length(seq)))
    cancer_muts_subset <- subset(cancer_muts_subset, aa_pos <= length(seq))
    cbind(cancer_muts_subset, data.frame(gene_symbol=gene_symbol, aln_pos=sapply(cancer_muts_subset$aa_pos, function(i) { which(i == gap_offsets) }), stringsAsFactors=F))
  }
  
  # FUNCTION 3 -----------------------------------------------------------------
  make_variant_matrix <- function(aln, variants, freq=F, threeletter=T) {
    variant_matrix <- matrix(0, nrow=length(aln), ncol=nchar(aln[1])) # make matrix of 0s
    rownames(variant_matrix) <- names(aln)
    
    # (1) LOOP #1 - MAKE EMPTY MATRIX FOR MASTER CHEMOKINE ALIGNMENT
    # acts on "aln" which corresponds to ccl_seqs
    for (s in rownames(variant_matrix)) {     # for each row...
      seq <- aln[s]                           # define new var that is a vector
                                              # seq positions from that row
      
      for (i in 1:(nchar(aln[1]))) {          # for each seq position
        if (substr(seq, i, i) == '-') {       # replace dashes (ie no AA in alignment)
          variant_matrix[s, i] <- NA          # with NA
        }
      }
    }                                         # after this loop, all matrix positions
                                              # with zeros indicate that a chemokine
                                              # has a residue there, and all NAs
                                              # indicate that a chemokine has a dash
    
    
    # (2) LOOP #2 - MULTIPLE MANIPULATIONS OF EMPTY (ZERO) MATRIX
    # acts on "aln" which corresponds to ccl_seqs AND acts on ccl_variants which
    # is the master SNP table from GNOMAD
    for(i in 1:nrow(variants)) {              # for each row in ccl_cancer_muts table
      r <- variants[i, ]                      # make vector of values from table
      seq <- aln[r$gene_symbol]               # grab gene name (eg CCL1) and seq
      
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
      
      # (2.3) IF-ELSE STATEMENT - IF "freq" IS TRUE
      if (freq) # if "freq" is TRUE...
        variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + r$MAF
      else
        variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + 1
    }
    
    variant_matrix
  }
  
  
##### 1: CANCER TABLE AND MATRIX GENERATION ####################################
  
  # (1.1) IMPORT FILES, SEQS FROM ENSEMBL --------------------------------------
  
    # import files, remove NAs, CXCR8
    cancer_muts <- read.csv("cancer/chemokinesMut5_20190207.csv", stringsAsFactors = F)
    cancer_muts <- cancer_muts %>% filter(!is.na(Mut_new))
    cancer_muts <- cancer_muts %>% filter(Gene != "GPR35") # remove CXCR8
    
    # fetch sequences based on Ensembl IDs from Daniel's table
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="http://feb2014.archive.ensembl.org")
    coding_seqs <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'ensembl_peptide_id', "ensembl_transcript_id", "coding"),
                         filters='ensembl_transcript_id', values=unique(cancer_muts$ENST_new), mart=ensembl)
    rm(ensembl)
    
    
  # (1.2) SELECT & ANNOTATE NONSENSE MUTATIONS WITHIN DANIEL'S TABLE -----------
    
    # select missense mutations and isolate NUCLEOTIE, FROM, and TO (2600 mutations)
    # Note that there are two syntaxes, so will treat individually and recombine
    cancer_muts <- cancer_muts %>% filter(Type == "Missense_Mutation")
    cancer_muts$temp_index <- c(1:nrow(cancer_muts))
    muts_parsed <- cancer_muts %>% 
      filter(grepl("c.[ACGT][0-9]+[ACGT]", Mut_new) |
               grepl("c.[0-9]+[ACGT]>[ACGT]", Mut_new)) %>%
      dplyr::select(Mut_new, temp_index)
    colnames(muts_parsed)[1] <- c("Mut")
    
    # isolate two different syntaxes, get NUCLEOTIDE, FROM, TO, & recombine
    synt.1 <- muts_parsed %>% filter(grepl("c.[ACGT][0-9]+[ACGT]", Mut)) # isolate
    synt.2 <- muts_parsed %>% filter(grepl("c.[0-9]+[ACGT]>[ACGT]", Mut)) # isolate
    synt.1.ann <- data.frame(str_match(synt.1$Mut, "c.([ACGT])([0-9]+)([ACGT])"), stringsAsFactors = F) # separate
    synt.2.ann <- data.frame(str_match(synt.2$Mut, "c.([0-9]+)([ACGT])>([ACGT])"), stringsAsFactors = F) # separate
    colnames(synt.1.ann) <- c("Mut_new","nt_from","pos", "nt_to") # annotate colnames
    colnames(synt.2.ann) <- c("Mut_new", "pos", "nt_from", "nt_to") # annotate colnames
    synt.1.ann$temp_index <- synt.1$temp_index # bring temporary index and NUC/FROM/TO together
    synt.2.ann$temp_index <- synt.2$temp_index # bring temporary index and NUC/FROM/TO together
    muts_parsed <- dplyr::bind_rows(synt.1.ann, synt.2.ann) # bring together fully annotated mutations from both syntaxes
    cancer_muts <- left_join(cancer_muts, muts_parsed, by = "temp_index") # bring together into original table
    cancer_muts$pos <- as.numeric(cancer_muts$pos)
    
    # remove used variables
    rm(muts_parsed, synt.1, synt.2, synt.1.ann, synt.2.ann)
    
  
  # (1.3) TRANSLATE NUCLEOTIDE TO AMINO ACID AND ADD TO DANIEL'S TABLE ---------
    
    # utilizes coding_seqs (SEQUENCE TABLE) and cancer_muts (mutation-annotated TCGA TABLE)
    cancer_muts_parsed <- adply(cancer_muts, 1, function(r) {
      print(r$ENST_new)
      coding_seq <- coding_seqs$coding[coding_seqs$ensembl_transcript_id == r$ENST_new]
      print(nchar(coding_seq))
      coding_split <- substring(coding_seq, seq(1, nchar(coding_seq)-2, 3), seq(3, nchar(coding_seq), 3))
      
      codon_coords <- idx_codon(r$pos)
      aa_pos <- codon_coords[1]
      codon_pos <- codon_coords[2]
      
      codon_from <- coding_split[aa_pos]
      codon_to <- codon_from
      substr(codon_to, codon_pos, codon_pos) <- r$nt_to
      
      print(paste(aa_pos, codon_from, codon_pos, r$nt_from))
      if(substr(codon_from, codon_pos, codon_pos) != r$nt_from) {
        print(paste("Mismatch in", r$gene_id, "at pos", aa_pos, codon_pos, "(expected", substr(codon_from, codon_pos, codon_pos), "saw", r$nt_from, ")"))
      }
      
      aa_ref <- GENETIC_CODE[codon_from]
      aa_to <- GENETIC_CODE[codon_to]
      
      data.frame(aa_pos=aa_pos, aa_ref=aa_ref, aa_to=aa_to, stringsAsFactors=F)
    })

    # remove used objects
    rm(cancer_muts, coding_seqs)
    
  
  # (1.4) MAP ANDY'S ALIGNMENT TO PARSED SEQUENCE POSITIONS --------------------
    
    # for dual-named chemokines/receptors, update table with conventional name
    cancer_muts_parsed <- cancer_muts_parsed %>% mutate(Gene = case_when(
      Gene == "PF4" ~ "CXCL4",
      Gene == "PF4V1" ~ "CXCL4L1",
      Gene == "PPBP" ~ "CXCL7",
      Gene == "IL8" ~ "CXCL8",
      Gene == "CX3CR1" ~ "CX3C1",
      Gene == "DARC" ~ "ACKR1",
      Gene != "PF4" | Gene != "PF4V1" | Gene != "PPBP" | Gene != "IL8" | Gene != "CX3CR1" | Gene != "DARC" ~ Gene
    ))
    
    # subset MUTATION TABLE for CHEMOKINES, add ALIGNMENT POSITION
    cancer_muts_parsed$gene_symbol <- cancer_muts_parsed$Gene # required for get_aligned_cancer_variants
    ccl_table <- read_csv("output/CK_SEQS.csv") # 
    ccl_cancer_muts <- adply(ccl_table, 1, get_aligned_cancer_variants) # match sequences
    
    # remove used objects
    rm(ccl_table)
    
    # import CHEMOKINE sequences
    aln.nterm <- read.fasta("sequences/chemokine_paralogs_NTERM.fasta", as.string=T)
    aln.core <- read.fasta("sequences/chemokine_paralogs_CORE.fasta", as.string=T)
    aln.cterm <- read.fasta("sequences/chemokine_paralogs_CTERM.fasta", as.string=T)
    aln.nterm <- aln.nterm[names(aln.core)]
    aln.cterm <- aln.cterm[names(aln.core)]
    ccl_seqs <- sapply(names(aln.core), function(n) { paste0(aln.nterm[[n]], aln.core[[n]], aln.cterm[[n]]) })
    gene_symbols <- toupper(str_split(names(aln.core), "_", simplify=T)[, 1])
    names(ccl_seqs) <- gene_symbols
    
    # import RECEPTOR sequences
    aln.ccr <- read.fasta("sequences/chemokine_receptor_paralogs.fasta", as.string=T)
    ccr_seqs <- sapply(aln.ccr, as.character)
    gene_symbols <- toupper(str_split(names(aln.ccr), "_", simplify=T)[, 1])
    names(ccr_seqs) <- gene_symbols
    ccr_gene_symbols <- toupper(str_split(names(aln.ccr), "_", simplify=T)[, 1])
    
    # subset MUTATION TABLE for CHEMOKINES, add ALIGNMENT POSITION
    ccr_table <- data.frame(gene_symbol=ccr_gene_symbols, aln=ccr_seqs, stringsAsFactors=F)
    ccr_cancer_muts <- adply(ccr_table, 1, get_aligned_cancer_variants)
    
    # remove used objects
    rm(aln.ccr, aln.core, aln.cterm, aln.nterm, ccr_table)
    
  # (1.4) CONVERT TO MATRIX & ADD NAMES ----------------------------------------
    
    # CHEMOKINE - make, write matrix
    cancer_muts_ccl_matrix <- make_variant_matrix(ccl_seqs, ccl_cancer_muts, threeletter=F)
    ccn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
    colnames(cancer_muts_ccl_matrix) <- c(ccn_names)
    cancer_muts_ccl_matrix <- cbind(rownames(cancer_muts_ccl_matrix), cancer_muts_ccl_matrix)
    colnames(cancer_muts_ccl_matrix)[1] <- c("protein")
    cancer_muts_ccl_matrix <- as.data.frame(cancer_muts_ccl_matrix)
    write_csv(cancer_muts_ccl_matrix, "output/CK_CANCER_MATRIX.csv")
    
    # RECEPTOR - make, write matrix
    cancer_muts_ccr_matrix <- make_variant_matrix(ccr_seqs, ccr_cancer_muts, threeletter=F)
    gn_names <- c(read.table("sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt", sep = ",", colClasses = "character"))
    colnames(cancer_muts_ccr_matrix) <- c(gn_names)
    cancer_muts_ccr_matrix <- cbind(rownames(cancer_muts_ccr_matrix), cancer_muts_ccr_matrix)
    colnames(cancer_muts_ccr_matrix)[1] <- c("protein")
    cancer_muts_ccr_matrix <- as.data.frame(cancer_muts_ccr_matrix)
    write_csv(cancer_muts_ccr_matrix, "output/CKR_CANCER_MATRIX.csv")
    
    # remove used objects
    rm(cancer_muts_ccl_matrix, cancer_muts_ccr_matrix, cancer_muts_parsed,
       ccl_cancer_muts, ccn_names, ccr_cancer_muts,
       gn_names, AMINO_ACID_CODE, ccl_seqs, ccr_seqs, gene_symbols, ccr_gene_symbols,
       get_aligned_cancer_variants, idx_codon, make_variant_matrix)
    
    