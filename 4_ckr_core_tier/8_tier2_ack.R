# Name:     2_tier2_chemokine_classification.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 7 (right panel)

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  library(ggridges)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/4_ckr_core_tier/")
  
##### FUNCTIONS ################################################################
  
  # FUNCTION 1 -----------------------------------------------------------------
  BuildTrainTest <- function(FILE, SEED, TRAIN.PCT){
    
    data <- read.csv(FILE) 
    data <- data %>% select(-protein, -seq)
    
    # define classes
    A <- data %>% filter(class == "cc") 
    B <- data %>% filter(class == "cxc")
    
    # random selection of training rows
    set.seed(SEED)  # for repeatability of samples
    rows.A <- sample(1:nrow(A), TRAIN.PCT*nrow(B))  # deliberate B for even sampling
    rows.B <- sample(1:nrow(B), TRAIN.PCT*nrow(B)) 
    train.A <- A[rows.A, ]  
    train.B <- B[rows.B, ]
    train <- rbind(train.A, train.B)  # row bind the 1's and 0's 
    
    # define test rows
    test.A <- A[-rows.A, ]
    test.B <- B[-rows.B, ] 
    test.A <- sample_n(test.A, nrow(test.B)) # match row numbers
    
    test <- rbind(test.A, test.B)  # row bind the 1's and 0's 
    train$cat <- c("train")
    test$cat <- c("test")
    all <- rbind(train, test)
    return(all)
    rm(A, B, rows.A, rows.B, train.A, train.B, all, test.A, test.B, all)
  }
  
  # FUNCTION 2 -----------------------------------------------------------------
  LogResAcc <- function(TRAIN, TEST){
    results <- data.frame()
    
    for(i in names(TRAIN)[2:ncol(TRAIN)]){
      
      fmla <- as.formula(paste0("class ~ ", i)) 
      logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
      logitMod$xlevels$i <- union(logitMod$xlevels$i, levels(test$i))
      noquote(paste0("logitMod$xlevels$", i, " <- union(logitMod$xlevels$", i, ", levels(test$", i, "))"))
      cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(logitMod$xlevels$", i, ", levels(test$", i, "))"))
      eval(parse(text = cmd))
      
      predicted <- plogis(predict(logitMod, TEST))  # predicted scores
      
      con.mat <- InformationValue::confusionMatrix(TEST$class, predicted) # see also ModelMetrics
      acc <- (con.mat[1,1] + con.mat[2,2]) / 
        (con.mat[1,1] + con.mat[2,2] + con.mat[1,2] + con.mat[2,1])
      acc <- as.data.frame(acc)
      acc$motif <- paste(i)
      results <- rbind(results, acc)
      rm(logitMod, predicted, con.mat, acc)
 
    }
    results <- results %>% arrange(desc(acc))
    return(results)
    rm(results)
  }
  
  # FUNCTION 3 -----------------------------------------------------------------
  LogResClass <- function(TRAIN){
    results <- data.frame()
    
    test.master <- read.csv("output/receptor_core_dataframe_ack_b2ar.csv") 
    test.master <- subset(test.master, grepl("human", test.master$seq))
    map.class <- test.master %>% select(protein, class)
    test.master$class <- as.character(test.master$class)
    protein.names <- test.master$protein
    protein.names <- as.character(protein.names)
    test.master$class <- as.character(test.master$class)
    
    test.master <- test.master %>% mutate(class = case_when(
      class == "ack" ~ "cc",
      class != "ack" ~ class
    )) %>% mutate(class = case_when(
      class == "cx3c" ~ "cc",
      class != "cx3c" ~ class
    )) %>% mutate(class = case_when(
      class == "xc" ~ "cc",
      class != "xc" ~ class
    )) %>% mutate(class = case_when(
      class == "non" ~ "cc",
      class != "non" ~ class
    ))
    
    for(j in protein.names){
      
      test <- test.master %>% filter(protein == j)
      test <- test %>% select(-protein, -seq)
      
      for(i in names(TRAIN)[2:ncol(TRAIN)]){
        
        fmla <- as.formula(paste0("class ~ ", i)) 
        logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
        logitMod$xlevels$i <- union(logitMod$xlevels$i, levels(test$i))
        noquote(paste0("logitMod$xlevels$", i, " <- union(logitMod$xlevels$", i, ", levels(test$", i, "))"))
        cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(logitMod$xlevels$", i, ", levels(test$", i, "))"))
        eval(parse(text = cmd))
        
        predicted <- plogis(predict(logitMod, test))  # predicted scores
        con.mat <- InformationValue::confusionMatrix(test$class, predicted) # see also ModelMetrics
        
        df <- as.data.frame(rownames(con.mat)[1])
        df$position <- paste(i)
        df$protein <- paste(j)
        results <- rbind(results, df)
        rm(logitMod, predicted, con.mat, df)
      }
      rm(test)
    }
    colnames(results)[1] <- c("cc_cxc")
    results$class <- map.class$class[match(unlist(results$protein), map.class$protein)]
    return(results)
    rm(results)
  }
  
  # FUNCTION 4 -----------------------------------------------------------------
  LogResScore <- function(TRAIN){
    results <- data.frame()
    
    test.master <- read.csv("output/receptor_core_dataframe_ack_b2ar.csv") 
    test.master <- subset(test.master, grepl("human", test.master$seq))
    map.class <- test.master %>% select(protein, class)
    test.master$class <- as.character(test.master$class)
    protein.names <- test.master$protein
    protein.names <- as.character(protein.names)
    test.master$class <- as.character(test.master$class)
    
    test.master <- test.master %>% mutate(class = case_when(
      class == "ack" ~ "cc",
      class != "ack" ~ class
    )) %>% mutate(class = case_when(
      class == "cx3c" ~ "cc",
      class != "cx3c" ~ class
    )) %>% mutate(class = case_when(
      class == "xc" ~ "cc",
      class != "xc" ~ class
    )) %>% mutate(class = case_when(
      class == "non" ~ "cc",
      class != "non" ~ class
    ))
    
    for(j in protein.names){
      
      test <- test.master %>% filter(protein == j)
      test <- test %>% select(-protein, -seq)
      
      for(i in names(TRAIN)[2:ncol(TRAIN)]){
        
        fmla <- as.formula(paste0("class ~ ", i)) 
        logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
        logitMod$xlevels$i <- union(logitMod$xlevels$i, levels(test$i))
        noquote(paste0("logitMod$xlevels$", i, " <- union(logitMod$xlevels$", i, ", levels(test$", i, "))"))
        cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(logitMod$xlevels$", i, ", levels(test$", i, "))"))
        eval(parse(text = cmd))
        
        predicted <- plogis(predict(logitMod, test))  # predicted scores
        con.mat <- InformationValue::confusionMatrix(test$class, predicted) # see also ModelMetrics
        
        df <- as.data.frame(predicted)
        df$position <- paste(i)
        df$protein <- paste(j)
        results <- rbind(results, df)
        rm(logitMod, predicted, con.mat, df)
      }
      rm(test)
    }
    colnames(results)[1] <- c("cc_cxc")
    results$class <- map.class$class[match(unlist(results$protein), map.class$protein)]
    return(results)
    rm(results)
  }
  
##### PART 1: TRAIN AND TEST WITH EACH RECEPTOR ################################

  # RUN 1
  # import, select training, test, run
  data <- BuildTrainTest("output/receptor_core_dataframe.csv", 100, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  accLR.1 <- LogResClass(train)
  accLR.1$run <- c("run1")
  
  scoreLR.1 <- LogResScore(train)
  scoreLR.1$run <- c("run1")
  
  # RUN 2
  # import, select training, test, run
  data <- BuildTrainTest("output/receptor_core_dataframe.csv", 200, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  accLR.2 <- LogResClass(train)
  accLR.2$run <- c("run2")
  
  scoreLR.2 <- LogResScore(train)
  scoreLR.2$run <- c("run2")
  
  # RUN 3
  # import, select training, test, run
  data <- BuildTrainTest("output/receptor_core_dataframe.csv", 200, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  accLR.3 <- LogResClass(train)
  accLR.3$run <- c("run3")
  
  scoreLR.3 <- LogResScore(train)
  scoreLR.3$run <- c("run3")
  
  # COMBINE
  acc <- bind_rows(accLR.1, accLR.2, accLR.3)
  write_csv(acc, "output/ackr_calc_acc.csv")
  scr <- bind_rows(scoreLR.1, scoreLR.2, scoreLR.3)
  write_csv(scr, "output/ackr_calc_score.csv")
  rm(accLR.1, accLR.2, accLR.3)
  rm(scoreLR.1, scoreLR.2, scoreLR.3)
  rm(train, data)
  
##### PART 2: ACCURACY AT CC & CXC CLASSIFICATION (BINARIZED OUTCOME) ##########
  
  # import and subset by Tier 2
  data <- read_csv("output/ackr_calc_acc.csv")
  tier2 <- read_csv("output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(mean >= 0.8) # %>% filter(int == "int")
  tier2 <- tier2$motif
  tier2 <- c(tier2, "gn3x32")
  data <- data %>% filter(position %in% tier2)
  data <- data %>% filter(position != "gnNTr.Cm10")
  
  # import interface
  rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
    filter(Chain1 != Chain2) %>% select(target_gnccn) %>% distinct()
  rin$gn <- c("gn")
  inter <- unite(rin, gn, sep = "", c(gn, target_gnccn))
  inter <- inter$gn
  inter <- c(inter,  "gn1x28_", "gn1x24_", "gn1x25_")
  rm(rin)

  data <- data %>% mutate(interface = case_when(
    position %in% inter ~ "yes",
    !(position %in% inter) ~ "no"
    ))
  data <- data %>% filter(interface == "yes")
  
  # score 0/1 by CC/CXC
  data <- data %>% mutate(cc_cxc = case_when(
    cc_cxc == 0 ~ "cc",
    cc_cxc == 1 ~ "cxc"
  ))
  
  # quantify 
  data <- data %>% 
    select(cc_cxc, protein, class, run) %>%
    count(cc_cxc, protein, class, run) %>%
    spread(cc_cxc, n)
  data[is.na(data)] <- 0
  
  data <- data %>%
    mutate(pct_cxc = cxc/14) # interface only
    #mutate(pct_cxc = cxc/39) # all T2
  
  data <- data %>% group_by(protein, class) %>% 
    mutate(mean = mean(pct_cxc)) %>% ungroup()
  
  # order and plot
  order <- as.factor(c("cc", "cxc", "cx3c", "xc", "ack", "non"))
  data$class <- factor(data$class, levels = order)
    
  data %>%
    filter(class != "xc") %>%
    filter(class != "cx3c") %>%
    ggplot(aes(class, pct_cxc)) +
    geom_violin() +
    #geom_point() +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, fill = "white", color="grey50") +
    stat_summary(fun.y=median, geom="point", shape=23, size=2) +
    theme_minimal()

  
##### PART 3: ACCURACY AT CC & CXC CLASSIFICATION BY SCORE #####################
  
  # import data
  data <- read_csv("output/ackr_calc_score.csv")
  
  # import and subset by Tier 2
  tier2 <- read_csv("output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(mean >= 0.8) # %>% filter(int == "int")
  tier2 <- tier2$motif
  tier2 <- c(tier2, "gn3x32")
  data <- data %>% filter(position %in% tier2)
  data <- data %>% filter(position != "gnNTr.Cm10")
  
  # import interface
  # rin <- read_csv("input/RIN.csv") %>% filter(class == "xray") %>%
  #   filter(Chain1 != Chain2) %>% select(target_gnccn) %>% distinct()
  # rin$gn <- c("gn")
  # inter <- unite(rin, gn, sep = "", c(gn, target_gnccn))
  # inter <- inter$gn
  # inter <- c(inter,  "gn1x28_", "gn1x24_", "gn1x25_")
  # rm(rin)
  # 
  # data <- data %>% mutate(interface = case_when(
  #   position %in% inter ~ "yes",
  #   !(position %in% inter) ~ "no"
  #   ))
  # data <- data %>% filter(interface == "yes")
  
  # mean value per receptor across three repeats
  mean.pos <- data %>% group_by(position, protein, class) %>% 
      # first get mean per pos (normal distribution for replicates)
    summarise(mean = mean(cc_cxc), sd = sd(cc_cxc)) %>% ungroup()
  colnames(mean.pos)[4] <- c("score")
  med.gpcr <- mean.pos %>% group_by(protein, class) %>% 
      # now get median per gpcr; non-normal
    summarise(median = median(score), sd = sd(score)) %>% ungroup()
  med.gpcr <- med.gpcr %>% filter(class != "cx3c") %>% filter(class != "xc")
  mean.pos <- mean.pos %>% filter(class != "cx3c") %>% filter(class != "xc")
  
  # order and plot by class
  order <- as.factor(c("cc", "cxc", "cx3c", "xc", "ack", "non"))
  med.gpcr$class <- factor(med.gpcr$class, levels = (order))
  mean.pos$class <- factor(mean.pos$class, levels = rev(order))
  
  med.gpcr %>%
    ggplot(aes(x = median, y = class)) +
    geom_density_ridges() +
    geom_vline(xintercept = 0.6) +
    geom_vline(xintercept = 0.4) +
    coord_flip() +
    theme_minimal()
  
  # mean.pos %>%
  #   ggplot(aes(x = score, y = class)) +
  #   geom_density_ridges() +
  #   geom_vline(xintercept = 0.6) +
  #   geom_vline(xintercept = 0.4) +
  #   theme_minimal()
  
  # mean value per receptor per position across three repeats
  mean.pos.all <- data %>% group_by(position, protein, class) %>% # first get mean per pos
    summarise(mean = mean(cc_cxc), sd = sd(cc_cxc)) %>% ungroup()
  
  proteins <- c("ccr1","ccr2", "ccr3", "ccr4",
                "ccr5","ccr6","ccr7","ccr8","ccr9",
                "ccr10","cxcr1","cxcr2","cxcr3",
                "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                "ackr1","ackr2","ackr3","ackr4",
                "ccrl2", "5ht1d", "nk1r","c5ar1")
  mean.pos.all <- mean.pos.all %>%
    filter(protein %in% proteins)
  
  # order and plot
  order <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                               "ccr5","ccr6","ccr7","ccr8","ccr9",
                               "ccr10","cxcr1","cxcr2","cxcr3",
                               "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                               "ackr1","ackr2","ackr3","ackr4",
                               "ccrl2", "5ht1d", "nk1r", "c5ar1")))  
  order <- as.factor(order)
  mean.pos.all$protein <- factor(mean.pos.all$protein, levels = (order))
  
  mean.pos.all %>%
    #filter(protein == "ackr1" | protein == "mc3r" | protein == "5ht1a") %>%
    filter(class != "xc" & class != "cx3c") %>%
    ggplot(aes(x = mean, y = protein)) +
    geom_density_ridges() +
    #geom_violin() +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=.1, fill = "white", color="grey50") +
    geom_hline(yintercept = 0.6) +
    geom_hline(yintercept = 0.4) +
    #geom_point(shape = 21,size = 2, stroke = 0.5, fill = "white") +
    #ylim(-0.5,1.5) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0)) +
    coord_flip()
  
  # per position by class
  # order and plot by class
  order <- as.factor(c("cc", "cxc", "cx3c", "xc", "ack", "non"))
  mean.pos.all$class <- factor(mean.pos.all$class, levels = (order))

  mean.pos.all %>%
    filter(class != "xc" & class != "cx3c") %>%
    ggplot(aes(x = mean, y = class)) +
    geom_density_ridges() +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0)) +
    coord_flip()
  
  
  
  