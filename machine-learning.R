library(tidyverse)
library(caret)
library(randomForest)

# Data sets: Rett, Dup15q, ASD
rettDmrFull <- read.delim("data/Individual/Rett_sig_individual_smoothed_DMR_methylation.txt", check.names = FALSE)
rettDmrFullCB <- read.delim("data/Consensus_background/Rett_consensus_background_individual_smoothed_DMR_methylation.txt", check.names = FALSE)
rettInfo <- read.csv("data/Sample_info/Rett_sample_info.csv")

dupDmrFull <- read.delim("data/Individual/Dup15q_sig_individual_smoothed_DMR_methylation.txt")
dupDmrFullCB <- read.delim("data/Consensus_background/Dup15_consensus_background_individual_smoothed_DMR_methylation.txt")
dupInfo <- read.csv("data/Sample_info/Dup15q_sample_info.csv")

asdDmrFull <- read.delim("data/Individual/ASD_sig_individual_smoothed_DMR_methylation.txt")
asdDmrFullCB <- read.delim("data/Consensus_background/ASD_consensus_background_individual_smoothed_DMR_methylation.txt")
asdInfo <- read.csv("data/Sample_info/ASD_sample_info.csv")

# Clean Dataset -------------------------------------------------------------
# exclude range of columns from 'width' to 'RawDiff', transpose
cleanData <- function(dmrFull) {
  data <- dmrFull %>% 
      as.tibble() %>% 
      select(-(width:RawDiff)) %>%
      unite(seqId1, seqnames, start, sep = ":") %>%
      unite(seqId, seqId1, end, sep = "-") %>%
      # transpose: cols to rows
      gather(sampleID, values, -seqId) %>% # cols to rows
      # transpose: rows to cols
      spread(seqId, values)
    return(data)
}

cleanDataCB <- function(dmrFull) {
  data <- dmrFull %>% 
    drop_na() %>%
    as.tibble() %>% 
    select(-(width:strand)) %>%
    unite(seqId1, seqnames, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") %>%
    # transpose: cols to rows
    gather(sampleID, values, -seqId) %>% # cols to rows
    # transpose: rows to cols
    spread(seqId, values)
  return(data)
}

# add diagnosis column, match and read in from file
# remove Sample ID column
cleanData2 <- function(dmrCleanData, sampleInfo) {
  dmrFinalData <- dmrCleanData %>% 
                  add_column(diagnosis = sampleInfo$Diagnosis[match(dmrCleanData$sampleID, sampleInfo$Name)], .after = 1) %>%
                  select(-sampleID)
  return(dmrFinalData)
}

rettDmr <- cleanData(rettDmrFull)
rettSampleID <- rettDmr$sampleID
rDmr <- cleanData2(rettDmr, rettInfo)
rettDmrCB <- cleanDataCB(rettDmrFullCB)
rDmrCB <- cleanData2(rettDmrCB, rettInfo)

dupDmr <- cleanData(dupDmrFull)
dupSampleID <- dupDmr$sampleID
dDmr <- cleanData2(dupDmr, dupInfo)
dupDmrCB <- cleanDataCB(dupDmrFullCB)
dDmrCB <- cleanData2(dupDmrCB, dupInfo)

asdDmr <- cleanData(asdDmrFull)
asdSampleID <- asdDmr$sampleID
aDmr <- cleanData2(asdDmr, asdInfo)
asdDmrCB <- cleanDataCB(asdDmrFullCB)
aDmrCB <- cleanData2(asdDmrCB, asdInfo)

# Partition data into training and testing --------------------------------
seed <- 9999
set.seed(seed)
trainIndex <- createDataPartition(dmrData$diagnosis, 
                                  p = 0.8,
                                  list = FALSE )

training <- dmrData[trainIndex, ]
testing <- dmrData[-trainIndex, ]


# Feature Selection -------------------------------------------------------
# https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# http://dataaspirant.com/2018/01/15/feature-selection-techniques-r/
# https://www.datacamp.com/community/tutorials/feature-selection-R-boruta


# Models ------------------------------------------------------------------
# Random forest, Neural networks, ant colony optimization
# particle swarm optimzation, genetic programming

# Model: Random Forest different trControl  ---------------------------
fitControl <- trainControl(method = "none", returnResamp = "final")
set.seed(seed)
rf_default <- train( diagnosis ~ ., 
                     data = training, 
                     method = "rf", 
                     trControl = fitControl)

# Model: Random Forest 10-fold Cross Validation ---------------------------
fitControl <- trainControl(method = "repeatedcv", # 10-fold cv
                           number = 10, 
                           repeats = 10, 
                           classProbs = TRUE) 
set.seed(seed)
rf_default <- train( diagnosis ~ ., 
                     data = training, 
                     method = "rf", 
                     trControl = fitControl)

# predict the outcome on a test set
rettPredict <- predict(rf_default, testing)
confusionMatrix(rettPredict, testing$diagnosis)

# Model: Random Forest using ranger ---------------------------------------
library(ranger)
rf_fit <- train(diagnosis ~ ., 
                data = training, 
                method = "ranger")
rf_fit

# predict the outcome on a test set
rangerPredict <- predict(rf_fit, testing)
confusionMatrix(rangerPredict, testing$diagnosis)

# Model: Stochastic Gradient Boosting -------------------------------------
library(gbm)
gbmFitControl <- trainControl(method = "repeatedcv", 
                              number = 2, 
                              repeats = 2)
set.seed(seed)
gbmFit1 <- train(diagnosis ~ ., 
                 data = training, 
                 method = "gbm", 
                 trControl = gbmFitControl, 
                 verbose = FALSE)

# Variable Importance for Random Forest -----------------------------------
vi <- varImp(object = rf_default)
vi
vi_plot <- plot(vi, main = "Random Forest - Variable Importance")
vi_plot

# warnings test: no warnings if no resampling -----------------------------------------------------------
# below fitControl and train give no resampling, no warnings
# but why no resampling? 
# warning ok? https://github.com/topepo/caret/issues/905

fitControl <- trainControl(method = "none", 
                           classProbs = TRUE, 
                           returnData = TRUE, 
                           returnResamp = "all", 
                           savePredictions = "all") 
fitControl

set.seed(seed)
rf_noresamp <- train( diagnosis ~ .,
                      data = training, 
                      method = "rf", 
                      trControl = fitControl) 
rf_noresamp 

noresampPredict <- predict(rf_noresamp, testing) 
confusionMatrix(noresampPredict, testing$diagnosis) 

# warnings: 
# 1. missing values in resampled performance measures
# caret/workflows.R
# resampleHist.R
# postResample.R
# resampleSummary.R
# resamples. R - collation and visualization of resampling results 

# caret:
# https://www.analyticsvidhya.com/blog/2016/12/practical-guide-to-implement-machine-learning-with-caret-package-in-r-with-practice-problem/