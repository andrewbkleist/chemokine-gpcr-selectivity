# Name:     2_tier2_chemokine_classification.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Supp. Figure 7 (right panel)

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/20190201_FINAL/4_ckr_core_tier/")
  
##### FUNCTIONS ################################################################
  
  # FUNCTION 1 -----------------------------------------------------------------
  BuildTrainTest <- function(FILE, SEED, TRAIN.PCT){
    
    data <- read.csv(FILE) 
    data <- data[, sapply(data, nlevels) > 1]
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
  
##### LOGISTIC REGRESSION ######################################################
  
  # RUN 1
  # import, select training, test, run
  data <- BuildTrainTest("output/receptor_core_dataframe.csv", 100, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  test <- data %>% filter(cat == "test") %>% select(-cat)
  accLR.ckr1 <- LogResAcc(train, test)
  rm(data, train, test)
  
  # RUN 2
  # import, split to training, test sets
  data <- BuildTrainTest("output/receptor_core_dataframe.csv", 200, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  test <- data %>% filter(cat == "test") %>% select(-cat)
  accLR.ckr2 <- LogResAcc(train, test)
  rm(data, train, test)
  
  # RUN 3
  # import, split to training, test sets
  data <- BuildTrainTest("output/receptor_core_dataframe.csv", 300, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  test <- data %>% filter(cat == "test") %>% select(-cat)
  accLR.ckr3 <- LogResAcc(train, test)
  rm(data, train, test)
  
  # COMBINE, STD DEV
  colnames(accLR.ckr1)[1] <- c("acc1")
  colnames(accLR.ckr2)[1] <- c("acc2")
  colnames(accLR.ckr3)[1] <- c("acc3")
  
  master <- left_join(accLR.ckr1, accLR.ckr2)
  master <- left_join(master, accLR.ckr3)
  master <- master %>% dplyr::select(motif, acc1, acc2, acc3)
  master <- master %>% gather(acc, value, 2:4)
  master <- master %>% group_by(motif) %>% summarise(mean = mean(value), sd = sd(value))
  master <- master %>% filter(!is.na(mean))
  
  # PLOT MEAN AND SD
  master$motif <- factor(master$motif, levels = master$motif[order(master$mean)])
  master %>%
    top_n(50, mean) %>%
    ggplot(aes(motif, mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()
  
  # write output
  write_csv(master, "output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
  