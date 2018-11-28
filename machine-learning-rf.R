library(tidyverse)
library(caret) 
library(randomForest)

# Datasets: Rett, Dup15q, ASD
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

info <- tibble(sampleID = c(as.character(rettInfo$Name), 
                            as.character(dupInfo$Name),
                            as.character(asdInfo$Name)), 
               diagnosis = c(as.character(rettInfo$Diagnosis), 
                             as.character(dupInfo$Diagnosis), 
                             as.character(asdInfo$Diagnosis)), 
               batch = c(as.character(rettInfo$batch), 
                         as.character(dupInfo$batch), 
                         as.character(asdInfo$batch)))

#' cleanData
#' @description Filter (exclude columns "width" to "RawDiff") and transpose DMR dataset
#' @param dmrFull DMR dataset 
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

#' cleanData2
#' @description Add diagnosis column for each sample by matching from sample info files and remove sample ID column
#' @param dmrCleanData dataset returned by cleanData()
#' @param sampleInfo sample info dataset
#' @import tidyverse
#' @export cleanData2
cleanData2 <- function(dmrCleanData, sampleInfo) {
  dmrFinalData <- dmrCleanData %>% 
    add_column(diagnosis = sampleInfo$Diagnosis[match(dmrCleanData$sampleID, sampleInfo$Name)], .after = 1) %>%
    select(-sampleID)
  return(dmrFinalData)
}

rettDmr <- cleanData(rettDmrFull)
rettSampleID <- rettDmr$sampleID
rDmr <- cleanData2(rettDmr, rettInfo)

dupDmr <- cleanData(dupDmrFull)
dupSampleID <- dupDmr$sampleID
dDmr <- cleanData2(dupDmr, dupInfo)

asdDmr <- cleanData(asdDmrFull)
asdSampleID <- asdDmr$sampleID
aDmr <- cleanData2(asdDmr, asdInfo)

#' cleanDataCB
#' @description Filter and transpose consensus background DMR dataset
#' @param dmrFull consensus background DMR dataset 
#' @import tidyverse
#' @export cleanDataCB
cleanDataCB <- function(dmrFull) {
  data <- dmrFull %>% 
    #drop_na() %>%
    as.tibble() %>% 
    select(-matches("width")) %>% 
    select(-matches("strand")) %>%
    unite(seqId1, seqnames, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") 
  return(data)
}

rettDmrCB <- cleanDataCB(rettDmrFullCB)
dupDmrCB <- cleanDataCB(dupDmrFullCB)
# remove repeated samples: JLKD063 = 1136 , JLKD066 = 1406, JLKD067 = 1711
asdDmrCB <- cleanDataCB(asdDmrFullCB) %>% select(-c("JLKD063", "JLKD066", "JLKD067"))

# joinedCB: combined CB data with diagnosis and batch
joinedCB <- rettDmrCB %>%
  full_join(dupDmrCB, by = "seqId") %>%
  full_join(asdDmrCB, by = "seqId") %>%
  drop_na() %>%
  gather(sampleID, values, -seqId) %>% # transpose: cols to rows
  spread(seqId, values) # transpose: rows to cols
joinedCB <- joinedCB %>%
  add_column(diagnosis = info$diagnosis[match(joinedCB$sampleID, info$sampleID)], .after = 1) %>%
  add_column(batch = info$batch[match(joinedCB$sampleID, info$sampleID)], .after = 2)
  


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
  return(list("confMat" = confMat, "preds" = predictModel))
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
  predConfMat <- predictConfMat(dmrPart, rfModel)
  result <- list("rfModel" = rfModel, "confMat" = predConfMat$confMat, "preds" = predConfMat$preds, "testingDiag" = dmrPart$testing$diagnosis)
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


# ROC curves --------------------------------------------------------------
# y axis: true positive rate (sensitivity)
# x axis: false positive rate (specificity)
# https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html
# https://stackoverflow.com/questions/30818188/roc-curve-in-r-using-rpart-package
library(ROCR)
# rDmrResult$confMat$byClass[1]

# confusionMatrix{caret} by default uses 0.5 cutoff. can't change?
# use confusion.matrix to change threshold values?
data(ROCR.simple)
pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels )
perf <- performance( pred, "tpr", "fpr" )
perf <- performance( pred, "sens", "spec" )
plot( perf )

library(rpart)
data("kyphosis", package = "rpart")
rp <- rpart(Kyphosis ~ ., data = kyphosis) #model
preds <- predict(rp, type = "prob")[,2]

# predictions are continuous predictions fo the classification
# labels are the binary truth for each variable
# pred <- prediction(as.vector(rDmrResult$preds), as.vector(rDmrResult$testingDiag))
pred <- prediction(c(1, 1), as.vector(rDmrResult$testingDiag))
perf <- performance(pred, "tpr", "fpr")
plot(perf, main = "ROC Curve for Random Forest")


# ROC Curve for rDMR
# use predict(rfModel, dmrPart$testing, type="prob")
# by default, cutoff = 0.5 and type = "raw", which gives a prediction: Control or Rett
# rf ROC curve example: http://scg.sdsu.edu/rf_r/
# https://machinelearningmastery.com/difference-test-validation-datasets/
# inputs to prediction(): predicted class probs for validation dataset
model <- aDmrResult$rfModel
dmrPart <- partitionData(rDmr)
testing <- dmrPart$testing

runFunctions <- function(dmrData) {
  dmrPart <- partitionData(dmrData)
  rfModel <- fitRandomForestModel(dmrPart$training)
  predConfMat <- predictConfMat(dmrPart, rfModel)
  result <- list("rfModel" = rfModel, "confMat" = predConfMat$confMat, "preds" = predConfMat$preds, "testingDiag" = dmrPart$testing$diagnosis)
  return(result)
}

predict(model, testing)
predictConfMat <- function(dmrPartData, fitModel) {
  predictModel <- predict(fitModel, dmrPartData$testing)
  confMat <- confusionMatrix(predictModel, dmrPartData$testing$diagnosis)
  return(list("confMat" = confMat, "preds" = predictModel))
}

# working 
dmrPart <- partitionData(rDmr)
rfModel <- fitRandomForestModel(dmrPart$training)
predict(rfModel, dmrPart$testing)
preds <- predict(rfModel, dmrPart$testing, type = "prob")
preds <- as.vector(preds[,2]) #Rett col
pred <- prediction(preds, as.vector(dmrPart$testing$diagnosis))
perf <- performance(pred, "tpr", "fpr")
plot(perf, main = "ROC Curve for Random Forest")



