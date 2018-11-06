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

# Models ------------------------------------------------------------------
# Random forest, Neural networks, ant colony optimization
# particle swarm optimzation, genetic programming

# Model: Random Forest different trControl  ---------------------------
fitControl <- trainControl(method = "none", returnResamp = "final")

fitControl <- trainControl(method = "repeatedcv", 
                           number = 2, 
                           repeats = 2, 
                           classProbs = TRUE) 

# Model: Random Forest 2-fold Cross Validation ---------------------------
set.seed(seed)
rf_model <- train( diagnosis ~ ., 
                   data = training, 
                   method = "rf", 
                   trControl = fitControl )
                  # preProcess = "nzv" makes no difference, resampling works

# predict the outcome on a test set
rfPredict <- predict(rf_model, testing)
confusionMatrix(rfPredict, testing$diagnosis)

# Model: Random Forest using ranger ---------------------------------------
library(ranger)
ranger_fit_model <- train(diagnosis ~ ., 
                data = training, 
                method = "ranger", 
                trControl = fitControl)
ranger_fit_model

# predict the outcome on a test set
rangerPredict <- predict(ranger_fit_model, testing)
confusionMatrix(rangerPredict, testing$diagnosis)

# Feature Selection -------------------------------------------------------
# https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# http://dataaspirant.com/2018/01/15/feature-selection-techniques-r/
# https://www.datacamp.com/community/tutorials/feature-selection-R-boruta

# FEATURE SELECTION - Variable Importance -----------------------------------
# after fitting model
set.seed(seed)
varImpList <- varImp(object = rf_model)
vi <- varImpList[[1]]
which(vi$Overall > 70)
vi[which(vi$Overall > 70), ]

dmrData_vi_rows <- row.names(vi)[which(vi$Overall > 70)]

vi_plot <- plot(vi, main = "Random Forest - Variable Importance")
vi_plot

# FEATURE SELECTION - Remove highly correlated variables
# before fitting model
# 0.75 cutoff: all models drastically worsen
# 0.85 cutoff: only asd better 
# 0.90 cutoff: all models slightly better
removeHighCor <- function(dmrData){
  set.seed(seed)
  cutoffValue = 0.90
  dmrData_noDiagnosis <- dmrData[, -1]
  corMatrix <- cor(dmrData_noDiagnosis)
  highCor <- findCorrelation(corMatrix, cutoff = cutoffValue)
  dmrData_noDiagnosis_noHC <- dmrData_noDiagnosis[, -highCor]
  dmrData_noHC <- add_column(dmrData_noDiagnosis_noHC, diagnosis = dmrData$diagnosis, .before = 1)
  return(dmrData_noHC)
}
rDmr_noHC <- removeHighCor(rDmr)
dDmr_noHC <- removeHighCor(dDmr)
aDmr_noHC <- removeHighCor(aDmr)

# FEATURE SELECTION - RFE recursive feature elimination
control <- rfeControl(functions = rfFuncs, 
                      method = "cv", 
                      number = 2)
# error: "need same number of samples in x and y" but they are the same
results <- rfe(x = training[,-1],
               y = training[, 1], 
               sizes = c(1:100), 
               rfeControl = control)
