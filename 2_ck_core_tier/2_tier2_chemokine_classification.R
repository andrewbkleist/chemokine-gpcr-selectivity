# Name:     2_chemokine_classification.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 2E

##### LOAD PACKAGES, SET WD ####################################################

  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  library(randomForest)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/2_ck_core_tier/")

  
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
  
  # note that some training/test set splits give error with some random seeds:
  # "Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
  # contrasts can be applied only to factors with 2 or more levels"
  # This is due to training/test set splits where there is no variation
  # in a specific residue; different random seeds were randomly selected to 
  # yield results with no errors on 3 runs and combined to get means and SD
  
  # RUN 1
  # import, select training, test, run
  # import, select training, test, run
  data <- BuildTrainTest("output/chemokine_core_dataframe.csv", 60, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  test <- data %>% filter(cat == "test") %>% select(-cat)
  accLR.ck1 <- LogResAcc(train, test)
  
  # RUN 2
  # import, split to training, test sets
  data <- BuildTrainTest("output/chemokine_core_dataframe.csv", 79, 0.8) #70
  train <- data %>% filter(cat == "train") %>% select(-cat)
  test <- data %>% filter(cat == "test") %>% select(-cat)
  accLR.ck2 <- LogResAcc(train, test)
  rm(data, train, test)
  
  # RUN 3
  # import, split to training, test sets
  data <- BuildTrainTest("output/chemokine_core_dataframe.csv", 101, 0.8)
  train <- data %>% filter(cat == "train") %>% select(-cat)
  test <- data %>% filter(cat == "test") %>% select(-cat)
  accLR.ck3 <- LogResAcc(train, test)
  rm(data, train, test)
  
  # COMBINE, STD DEV
  colnames(accLR.ck1)[1] <- c("acc1")
  colnames(accLR.ck2)[1] <- c("acc2")
  colnames(accLR.ck3)[1] <- c("acc3")
  
  master <- left_join(accLR.ck1, accLR.ck2)
  master <- left_join(master, accLR.ck3)
  master <- master %>% dplyr::select(motif, acc1, acc2, acc3)
  master <- master %>% gather(acc, value, 2:4)
  master <- master %>% group_by(motif) %>% summarise(mean = mean(value), sd = sd(value))
  master <- master %>% filter(!is.na(mean))
  
  # WRITE OUTPUT
  write_csv(master, "output/CHEMOKINE_CLASSIFICATION_N3.csv")
  
  # PLOT MEAN AND SD
  master$motif <- factor(master$motif, levels = master$motif[order(master$mean)])
  master %>%
    #top_n(20, mean) %>%
    filter(mean >= .8) %>%
    ggplot(aes(motif, mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()
  
  