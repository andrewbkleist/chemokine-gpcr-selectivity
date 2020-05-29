# Name:     1_generate_snp.R
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
  AMINO_ACID_CODE <- toupper(AMINO_ACID_CODE)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/6_snp_cancer/")
  
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
  
  
##### 1: SNP TABLE AND MATRIX GENERATION #######################################
#
#   SUMMARY OF DATA GENERATED
#
#   ccl_table               = sequence string for each chemokine with chemokine name
#   ccl_variant_matrix      = for each chemokine (rows) by each positions (cols)
#                             gives number of variants for that chemokine/position
#   ccl_variant_matrix      = same as above for chemokine core
#   ccl_variant_matrix_freq = same as above *matrix except freq not counts
#   ccl_variants            = variants from Gnomad, including reference, alternate,
#                             chromosome, position, annotation, allele no., allele
#                             count, allele frequency, etc.
#
#   same data frames are generated for receptors
  
  
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
    
    # KEY STEP 1 - FUNCTION 1 - GNOMAD to TABLE
    ccl_variants <- adply(ccl_table, 1, get_aligned_variants) 
      
    # KEY STEP 2 - FUNCTION 2 - GNOMAD to ALN MATRIX (no. muts at each pos)
    ccl_variant_matrix <- make_variant_matrix(ccl_seqs, ccl_variants)
    ccl_variant_matrix_freq <- make_variant_matrix(ccl_seqs, ccl_variants, allele.count =T)

    # removed used objects
    rm(gene_symbols, ccl_seqs, aln.core) # chemokine objects
    
    
  # (1.3) RECEPTOR - PARSE GNOMAD - RAW DATA ----------------------------------
    
    # get receptor GENE NAMES
    ccr_gene_symbols <- toupper(str_split(names(aln.ccr), "_", simplify=T)[, 1])
    names(ccr_seqs) <- ccr_gene_symbols
    ccr_table <- data.frame(gene_symbol=ccr_gene_symbols, aln=ccr_seqs, stringsAsFactors=F)
    
    # KEY STEP 1 - FUNCTION 1 - GNOMAD to TABLE
    ccr_variants <- adply(ccr_table, 1, get_aligned_variants)

    # KEY STEP 2 - FUNCTION 2 - GNOMAD to ALN MATRIX (no. muts at each pos)
    ccr_variant_matrix <- make_variant_matrix(ccr_seqs, ccr_variants)
    ccr_variant_matrix_freq <- make_variant_matrix(ccr_seqs, ccr_variants, allele.count =T)

    # removed used objects
    rm(aln.ccr, ccr_gene_symbols, ccr_seqs) # receptor objects
    rm(AMINO_ACID_CODE, format_variants_gnomad, get_aligned_variants, make_variant_matrix)
      
    
  # (1.3) TIDY AND WRITE OUTPUT ------------------------------------------------
    
    # convert to data frame
    ccl_variant_matrix <- as.data.frame(ccl_variant_matrix)
    ccl_variant_matrix_freq <- as.data.frame(ccl_variant_matrix_freq)
    ccr_variant_matrix <- as.data.frame(ccr_variant_matrix)
    ccr_variant_matrix_freq <- as.data.frame(ccr_variant_matrix_freq)
    
    # add names
    ccn_names <- c(read.table("sequences/FULL_CHEMOKINE_CCN.txt", sep = ",", colClasses = "character"))
    colnames(ccl_variant_matrix) <- c(ccn_names)
    colnames(ccl_variant_matrix_freq) <- c(ccn_names)
    gn_names <- c(read.table("sequences/FULL_RECEPTOR_GN_UNIQUE_GN.txt", sep = ",", colClasses = "character"))
    colnames(ccr_variant_matrix) <- c(gn_names)
    colnames(ccr_variant_matrix_freq) <- c(gn_names)
    
    # remove objects
    rm(gn_names, ccn_names)
    
    # add chemokine/receptor names
    ccl_variant_matrix <- cbind(rownames(ccl_variant_matrix), ccl_variant_matrix)
    colnames(ccl_variant_matrix)[1] <- c("protein")
    ccr_variant_matrix <- cbind(rownames(ccr_variant_matrix), ccr_variant_matrix)
    colnames(ccr_variant_matrix)[1] <- c("protein")
    ccl_variant_matrix_freq <- cbind(rownames(ccl_variant_matrix_freq), ccl_variant_matrix_freq)
    colnames(ccl_variant_matrix_freq)[1] <- c("protein")
    ccr_variant_matrix_freq <- cbind(rownames(ccr_variant_matrix_freq), ccr_variant_matrix_freq)
    colnames(ccr_variant_matrix_freq)[1] <- c("protein")
    
    # write output
    write_csv(ccl_table, "output/CK_SEQS.csv")
    write_csv(ccl_variants, "output/CK_SNP_GNOMAD.csv")
    write_csv(ccl_variant_matrix, "output/CK_SNP_MATRIX.csv")
    write_csv(ccl_variant_matrix_freq, "output/CK_SNP_FREQ_MATRIX.csv")
    
    write_csv(ccr_table, "output/CKR_SEQS.csv")
    write_csv(ccr_variants, "output/CKR_SNP_GNOMAD.csv")
    write_csv(ccr_variant_matrix, "output/CKR_SNP_MATRIX.csv")
    write_csv(ccr_variant_matrix_freq, "output/CKR_SNP_FREQ_MATRIX.csv")

    rm(ccl_variant_matrix, ccl_variant_matrix_freq, ccl_variants, ccl_table,
       ccr_variant_matrix, ccr_variant_matrix_freq, ccr_variants, ccr_table)

    