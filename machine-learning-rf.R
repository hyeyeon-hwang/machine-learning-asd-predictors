library(tidyverse)
library(caret) 
library(randomForest)

# Data sets: Rett, Dup15q, ASD
rettDmrFull <- read.delim("../data/Individual/Rett_sig_individual_smoothed_DMR_methylation.txt", check.names = FALSE)
rettDmrFullCB <- read.delim("../data/Consensus_background/Rett_consensus_background_individual_smoothed_DMR_methylation.txt", check.names = FALSE)
rettInfo <- read.csv("../data/Sample_info/Rett_sample_info.csv") 
rettInfo <- rettInfo %>% add_column(batch = "batchRett")

dupDmrFull <- read.delim("../data/Individual/Dup15q_sig_individual_smoothed_DMR_methylation.txt")
dupDmrFullCB <- read.delim("../data/Consensus_background/Dup15_consensus_background_individual_smoothed_DMR_methylation.txt")
dupInfo <- read.csv("../data/Sample_info/Dup15q_sample_info.csv") 
dupInfo <- dupInfo %>% add_column(batch = "batchDup")

asdDmrFull <- read.delim("../data/Individual/ASD_sig_individual_smoothed_DMR_methylation.txt")
asdDmrFullCB <- read.delim("../data/Consensus_background/ASD_consensus_background_individual_smoothed_DMR_methylation.txt")
asdInfo <- read.csv("../data/Sample_info/ASD_sample_info.csv")
asdInfo <- asdInfo %>% add_column(batch = "batchAsd")

sampleInfo <- tibble(sampleID = c(as.character(rettInfo$Name), 
                            as.character(dupInfo$Name),
                            as.character(asdInfo$Name)), 
               diagnosis = c(as.character(rettInfo$Diagnosis), 
                             as.character(dupInfo$Diagnosis), 
                             as.character(asdInfo$Diagnosis)), 
               batch = c(as.character(rettInfo$batch), 
                         as.character(dupInfo$batch), 
                         as.character(asdInfo$batch)))

#' cleanData
#' @description Filter (exclude columns "width" to "RawDiff") and transpose dmr dataset
#' @param Dmr dataset 
#' @import tidyverse
#' @export cleanData
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
    #drop_na() %>%
    as.tibble() %>% 
    select(-matches("width")) %>% 
    select(-matches("strand")) %>%
    unite(seqId1, seqnames, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") 
    # transpose: cols to rows
    ##gather(sampleID, values, -seqId) %>% # cols to rows
    # transpose: rows to cols
    ##spread(seqId, values)
  return(data)
}

# add diagnosis column, match and read in from file
# remove Sample ID column
cleanData2 <- function(dmrCleanData, sampleInfo) {
  dmrFinalData <- dmrCleanData %>% 
    add_column(diagnosis = sampleInfo$Diagnosis[match(dmrCleanData$sampleID, sampleInfo$Name)], .after = 1) %>%
    add_column(batch = )
    select(-sampleID)
  return(dmrFinalData)
}

rettDmrCB <- cleanDataCB(rettDmrFullCB)
dupDmrCB <- cleanDataCB(dupDmrFullCB)
# repeated samples: JLKD063 = 1136 , JLKD066 = 1406, JLKD067 = 1711
asdDmrCB <- cleanDataCB(asdDmrFullCB) %>% select(-c("JLKD063", "JLKD066", "JLKD067"))

# combined CB data with diagnosis and batch
joinedCB <- rettDmrCB %>%
  full_join(dupDmrCB, by = "seqId") %>%
  full_join(asdDmrCB, by = "seqId") %>%
  drop_na() %>%
  gather(sampleID, values, -seqId) %>% # cols to rows
  spread(seqId, values) %>% #rows to cols
  add_column(diagnosis = sampleInfo$diagnosis[match(joinedCB$sampleID, sampleInfo$sampleID)], .after = 1) %>%
  add_column(batch = sampleInfo$batch[match(joinedCB$sampleID, sampleInfo$sampleID)], .after = 2)
  
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



seed <- 9999
# Partition data into training and testing --------------------------------
partitionData <- function(dmrDataIn) {
  set.seed(seed)
  trainIndex <- createDataPartition(dmrDataIn$diagnosis, 
                                    p = 0.8,
                                    list = FALSE )
  
  training <- dmrDataIn[trainIndex, ]
  testing <- dmrDataIn[-trainIndex, ]
  dmrDataOut <- list("training" = training, "testing" = testing)
  return(dmrDataOut)
}

# Models ------------------------------------------------------------------
# Random forest, Neural networks, ant colony optimization
# particle swarm optimzation, genetic programming

# Model: Random Forest different trControl  ---------------------------
fitRandomForestModel <- function(trainingData) {
  #fitControl <- trainControl(method = "none", returnResamp = "final")
  fitControl <- trainControl(method = "repeatedcv", 
                             number = 2, 
                             repeats = 2, 
                             classProbs = TRUE) 
  
  # Model: Random Forest 2-fold Cross Validation ---------------------------
  set.seed(seed)
  rf_model <- train( diagnosis ~ ., 
                     data = trainingData, 
                     method = "rf", 
                     trControl = fitControl )
  # preProcess = "nzv" makes no difference, resampling works
  
  return(rf_model)
}

# predict the outcome on a test set
predictConfMat <- function(dmrPartData, fitModel) {
  predictModel <- predict(fitModel, dmrPartData$testing)
  confMat <- confusionMatrix(predictModel, dmrPartData$testing$diagnosis)
  return(confMat)
}

# Feature Selection -------------------------------------------------------
# https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# http://dataaspirant.com/2018/01/15/feature-selection-techniques-r/
# https://www.datacamp.com/community/tutorials/feature-selection-R-boruta

# FEATURE SELECTION - Remove highly correlated variables
# before fitting model
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

# FEATURE SELECTION - Variable Importance -----------------------------------
# after fitting model
selectImpVar <- function(dmrData, rf_model, cutoffValue = 70) {
  set.seed(seed)
  varImpList <- varImp(object = rf_model)
  vi <- varImpList[[1]]
  dmrData_vi_rows <- row.names(vi)[which(vi$Overall > cutoffValue)] 
  dmrData_vi_rows <- gsub("`", "", dmrData_vi_rows)
  dmrData_vi <- dmrData[, dmrData_vi_rows]
  dmrData_vi <- add_column(dmrData_vi, diagnosis = dmrData$diagnosis, .before = 1)
  return(dmrData_vi)
  #vi_plot <- plot(dmrData_vi, main = "Random Forest - Variable Importance")
  #vi_plot
}

# FEATURE SELECTION - RFE recursive feature elimination
control <- rfeControl(functions = rfFuncs, 
                      method = "cv", 
                      number = 2)
# error: "need same number of samples in x and y" but they are the same
results <- rfe(x = training[,-1],
               y = training[, 1], 
               sizes = c(1:100), 
               rfeControl = control)

# Run ---------------------------------------------------------------------

runFunctions <- function(dmrData) {
  dmrPart <- partitionData(dmrData)
  rfModel <- fitRandomForestModel(dmrPart$training)
  confMat <- predictConfMat(dmrPart, rfModel)
  result <- list("rfModel" = rfModel, "confMat" = confMat)
  return(result)
}

rDmrResult <- runFunctions(rDmr)
dDmrResult <- runFunctions(dDmr)
aDmrResult <- runFunctions(aDmr)

# run after selecting important variables
rDmr_vi <- selectImpVar(rDmr, rDmrResult$rfModel, cutoffValue = 70) # accuracy: 1 -> 0.5
dDmr_vi <- selectImpVar(dDmr, dDmrResult$rfModel, cutoffValue = 70) # accuracy: 1 -> 1
aDmr_vi <- selectImpVar(aDmr, aDmrResult$rfModel, cutoffValue = 70) # accuracy: 0.8 -> 1

rDmr_vi <- selectImpVar(rDmr, rDmrResult$rfModel, cutoffValue = 80) # accuracy: 1 -> 0.5
dDmr_vi <- selectImpVar(dDmr, dDmrResult$rfModel, cutoffValue = 80) # accuracy: 1 -> 1
aDmr_vi <- selectImpVar(aDmr, aDmrResult$rfModel, cutoffValue = 80) # accuracy: 0.8 -> 1

rDmrResult_vi <- runFunctions(rDmr_vi) 
dDmrResult_vi <- runFunctions(dDmr_vi) 
aDmrResult_vi <- runFunctions(aDmr_vi) 

# run after removing highly correlated variables
rDmr_noHC <- removeHighCor(rDmr)
dDmr_noHC <- removeHighCor(dDmr)
aDmr_noHC <- removeHighCor(aDmr)

# 0.75 cutoff: all models drastically worsen
# 0.85 cutoff: only asd better 
# 0.90 cutoff: all models slightly better
rDmrResult_noHC <- runFunctions(rDmr_noHC)
dDmrResult_noHC <- runFunctions(dDmr_noHC)
aDmrResult_noHC <- runFunctions(aDmr_noHC)
