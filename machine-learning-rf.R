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
dmrData <- rDmr
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
rf_model <- train( diagnosis ~ ., 
                     data = training, 
                     method = "rf", 
                     trControl = fitControl )

# Model: Random Forest 10-fold Cross Validation ---------------------------
fitControl <- trainControl(method = "repeatedcv", # 10-fold cv
                           number = 2, 
                           repeats = 2, 
                           classProbs = TRUE) 
set.seed(seed)

# removing near zero variance variables 
rf_model <- train( diagnosis ~ ., 
                     data = training, 
                     method = "rf", 
                     trControl = fitControl, 
                     preProcess = "nzv" ) #resampling works with or without

# predict the outcome on a test set
rettPredict <- predict(rf_model, testing)
confusionMatrix(rettPredict, testing$diagnosis)

# Model: Random Forest using ranger ---------------------------------------
library(ranger)
rf_fit <- train(diagnosis ~ ., 
                data = training, 
                method = "ranger", 
                trControl = fitControl)
rf_fit

# predict the outcome on a test set
rangerPredict <- predict(rf_fit, testing)
confusionMatrix(rangerPredict, testing$diagnosis)

# FEATURE SELECTION - Variable Importance -----------------------------------
set.seed(seed)
vi <- varImp(object = rf_model)
vi
vi_plot <- plot(vi, main = "Random Forest - Variable Importance")
vi_plot

# remove highly correlated variables
# 4641 variables total, 4502 highly correlated
set.seed(seed)
correlationMatrix <- cor(rDmr[,-1])
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.75)

# RFE: recursive feature elimination
control <- rfeControl(functions = rfFuncs, 
                      method = "cv", 
                      number = 2)

# error: need same number of samples in x and y, but they are the same
results <- rfe(x = training[,-1],
               y = training[, 1], 
               sizes = c(1:100), 
               rfeControl = control)

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