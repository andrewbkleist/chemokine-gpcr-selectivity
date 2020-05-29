# Name:     4_motif_log_reg.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 5B (right)

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  library(ggrepel)
  library(InformationValue)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/3_ck_motif/")

##### FUNCTIONS ################################################################
  
  # FUNCTION 1 -----------------------------------------------------------------
  BuildTrainTest <- function(FILE, SEED, TRAIN.PCT){
    data <- read_csv(FILE) 
      # Note that "NA" (Asn-Ala) from 2mer list is read as column name X224 or X3928
    
    # define classes
    A <- data %>% filter(class == "cc") 
    B <- data %>% filter(class == "cxc")
    
    # random selection of training rows
    set.seed(SEED)  # for repeatability of samples
    rows.A <- sample(1:nrow(A), TRAIN.PCT*nrow(B))  # cxc is deliberate here to make sampled rows equal
    rows.B <- sample(1:nrow(B), TRAIN.PCT*nrow(B)) 
    train.A <- A[rows.A, ]  
    train.B <- B[rows.B, ]
    train <- rbind(train.A, train.B)  # row bind the 1's and 0's 
    
    # define test rows
    test.A <- A[-rows.A, ]
    test.B <- B[-rows.B, ] 
    test.A <- sample_n(test.A, nrow(test.B)) # match row numbers
    
    test <- rbind(test.A, test.B)  # row bind the 1's and 0's 
    return(list(train, test))
    rm(A, B, rows.A, rows.B, train.A, train.B, all, test.A, test.B)
  }
  
  
  # FUNCTION 2 -----------------------------------------------------------------
  AccuracyEval <- function(TRAIN, TEST){
    results <- data.frame()
    TRAIN <- as.data.frame(TRAIN)
    TEST <- as.data.frame(TEST)
    
    TRAIN[,1] <- as.factor(TRAIN[,1])
    TEST[,1] <- as.factor(TEST[,1])
    
    for(i in names(TRAIN)[2:ncol(TRAIN)]){
      
      fmla <- as.formula(paste0("class ~ ", i)) 
      logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
      predicted <- plogis(predict(logitMod, TEST))  # predicted scores
      
      #optCutOff <- optimalCutoff(test.data$class, predicted)[1] 
      con.mat <- confusionMatrix(TEST$class, predicted)
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
  
##### 1: TRAIN MODEL WITH FEATURES INDIVIDUALLY ################################
  
  # RUN 1
  # import, split to training, test sets
  mer <- BuildTrainTest("output/CK_MOTIF_MATRIX.csv", 100, 0.8)
  train <- mer[[1]]
  test <- mer[[2]]
  train <- train %>% select(-file, -protein)
  test <- test %>% select(-file, -protein)
  
  # evaluate
  accuracy1 <- AccuracyEval(train, test)
  rm(mer, train, test)
  
  # RUN 2
  # import, split to training, test sets
  mer <- BuildTrainTest("output/CK_MOTIF_MATRIX.csv", 200, 0.8)
  train <- mer[[1]]
  test <- mer[[2]]
  train <- train %>% select(-file, -protein)
  test <- test %>% select(-file, -protein)
  
  # evaluate
  accuracy2 <- AccuracyEval(train, test)
  rm(mer, train, test)
  
  # RUN 3
  # import, split to training, test sets
  mer <- BuildTrainTest("output/CK_MOTIF_MATRIX.csv", 300, 0.8)
  train <- mer[[1]]
  test <- mer[[2]]
  train <- train %>% select(-file, -protein)
  test <- test %>% select(-file, -protein)
  
  # evaluate
  accuracy3 <- AccuracyEval(train, test)
  rm(mer, train, test)
  
  # COMBINE, STD DEV
  colnames(accuracy1)[1] <- c("acc1")
  colnames(accuracy2)[1] <- c("acc2")
  colnames(accuracy3)[1] <- c("acc3")
  
  master <- left_join(accuracy1, accuracy2)
  master <- left_join(master, accuracy3)
  master <- master %>% dplyr::select(motif, acc1, acc2, acc3)
  master <- master %>% gather(acc, value, 2:4)
  master <- master %>% group_by(motif) %>% summarise(mean = mean(value), sd = sd(value))
  master <- master %>% filter(!is.na(mean))
  
  # PLOT MEAN AND SD
  master$motif <- factor(master$motif, levels = master$motif[order(master$mean)])
  master %>%
    top_n(20, mean) %>%
    ggplot(aes(motif, mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()
  