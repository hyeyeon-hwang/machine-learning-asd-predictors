source("asd-predictors.R")
library(knitr)
library(kableExtra)

# Random Forest Model Results ---------------------------------------------
rDmrResult <- runFunctions(rDmr, p = 0.8, pos = "Rett")
dDmrResult <- runFunctions(dDmr, p = 0.8, pos = "Dup15q")
aDmrResult <- runFunctions(aDmr, p = 0.8, pos = "ASD")
pDmrResult <- runFunctions(pDmr, p = 0.8, pos = "idiopathic_autism")
# pDmrCBResult <- runFunctions(pDmrCB)

# confusion matrix tables using kabble
cmTable <- function(dmrResult) {
  dmrResult$confMat$table %>%
    kable() %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), font_size = 15) %>%
    add_header_above(header = c("Confusion Matrix" = 3), align = "c")
}
cmTable(rDmrResult)
cmTable(dDmrResult)
cmTable(aDmrResult)
cmTable(pDmrResult)

sumRes <- function(dmrResult, caption) {
  a <- paste("Model", dmrResult$rfModel$modelInfo$label, sep = ": ")
  b <- paste("Disorder", dmrResult$confMat$positive, sep = ": ")
  c <- paste(a, b, sep = ", ")
  ntrain <- dim(dmrResult$rfModel$trainingData)[1]
  ntest <- nrow(dmrResult$probPreds)
  npred <- dim(dmrResult$rfModel$trainingData)[2] - 1
  
  sumTable <- tibble(Measure = as.character(), Value = as.numeric()) %>%
    add_row(Measure = "Number of Samples in Training Data", Value = round(ntrain)) %>%
    add_row(Measure = "Number of Samples in Testing Data", Value = ntest) %>%
    add_row(Measure = "Number of Predictors", Value = npred) %>%
    add_row(Measure = "Accuracy", Value = dmrResult$confMat$overall["Accuracy"]) %>%
    add_row(Measure = "Kappa", Value = dmrResult$confMat$overall["Kappa"]) %>%
    add_row(Measure = "Accuracy P Value (Acc > NIR)", Value = dmrResult$confMat$overall["AccuracyPValue"]) %>%
    add_row(Measure = "Sensitivity", Value = dmrResult$confMat$byClass["Sensitivity"]) %>%
    add_row(Measure = "Specificity", Value = dmrResult$confMat$byClass["Specificity"]) %>%
    add_row(Measure = "Positive Predictive Values", Value = dmrResult$confMat$byClass["Pos Pred Value"]) %>%
    add_row(Measure = "Negative Predictive Values", Value = dmrResult$confMat$byClass["Neg Pred Value"]) %>%
    kable(caption = paste(c, caption, sep = " - ")) %>%
    kable_styling() %>%
    column_spec(1:2, color = "black") %>%
    add_header_above(header = c("Summarized results from classification algorithm" = 2), 
                     align = "c")
  return(sumTable)
}

rf_rett_table <- sumRes(rDmrResult, caption = "")
rf_rett_table
rf_dup_table <- sumRes(dDmrResult, caption = "")
rf_dup_table
rf_asd_table <- sumRes(aDmrResult, caption = "")
rf_asd_table
rf_plac_table <- sumRes(pDmrResult, caption = "")
rf_plac_table

# Feature Selection - Remove highly correlated variables ------------------
# run after removing highly correlated variables
rDmr_noHC_90 <- removeHighCor(rDmr, cutoffValue = 0.90) #10 samples, 2884 predictors, acc = 0.925, 1
rDmr_noHC_80 <- removeHighCor(rDmr, cutoffValue = 0.80) #10 samples, 460 predictors, acc = 0.2027, 0
rDmr_noHC_70 <- removeHighCor(rDmr, cutoffValue = 0.70)

dDmr_noHC_90 <- removeHighCor(dDmr, cutoffValue = 0.90) #8 samples, 1196 predictors, acc = 0.769, 0.5
dDmr_noHC_80 <- removeHighCor(dDmr, cutoffValue = 0.80) #8 samples, 116 predictors, acc = 0.0694, 0
dDmr_noHC_70 <- removeHighCor(dDmr, cutoffValue = 0.70)

aDmr_noHC_90 <- removeHighCor(aDmr, cutoffValue = 0.90) #22 samples, 470 predictors, acc = 0.849, 1
aDmr_noHC_80 <- removeHighCor(aDmr, cutoffValue = 0.80) #22 samples, 443 predictors, acc = 0.87, 1
aDmr_noHC_70 <- removeHighCor(aDmr, cutoffValue = 0.70)

pDmr_noHC_90 <- removeHighCor(pDmr, cutoffValue = 0.90) 
pDmr_noHC_80 <- removeHighCor(pDmr, cutoffValue = 0.80) 
pDmr_noHC_70 <- removeHighCor(pDmr, cutoffValue = 0.70)
pDmr_noHC_75 <- removeHighCor(pDmr, cutoffValue = 0.75)
pDmr_noHC_79 <- removeHighCor(pDmr, cutoffValue = 0.79)
pDmr_noHC_60 <- removeHighCor(pDmr, cutoffValue = 0.60)

# random forest models on data after removing highly correlated variables
rDmrResult_noHC_90 <- runFunctions(rDmr_noHC_90, p = 0.8, pos = "Rett") 
rDmrResult_noHC_80 <- runFunctions(rDmr_noHC_80, p = 0.8, pos = "Rett") 
rDmrResult_noHC_70 <- runFunctions(rDmr_noHC_70, p = 0.8, pos = "Rett")

dDmrResult_noHC_90 <- runFunctions(dDmr_noHC_90, p = 0.8, pos = "Dup15q")
dDmrResult_noHC_80 <- runFunctions(dDmr_noHC_80, p = 0.8, pos = "Dup15q")
dDmrResult_noHC_70 <- runFunctions(dDmr_noHC_70, p = 0.8, pos = "Dup15q")

aDmrResult_noHC_90 <- runFunctions(aDmr_noHC_90, p = 0.8, pos = "ASD")
aDmrResult_noHC_80 <- runFunctions(aDmr_noHC_80, p = 0.8, pos = "ASD")
aDmrResult_noHC_70 <- runFunctions(aDmr_noHC_70, p = 0.8, pos = "ASD")

# error, no predictors for 0.90, 0.80 cutoffs
pDmrResult_noHC_90 <- runFunctions(pDmr_noHC_90, p = 0.8, pos = "idiopathic_autism") #no predictors
pDmrResult_noHC_80 <- runFunctions(pDmr_noHC_80, p = 0.8, pos = "idiopathic_autism") #no predictors
pDmrResult_noHC_70 <- runFunctions(pDmr_noHC_70, p = 0.8, pos = "idiopathic_autism") # less 1, 286 predictors
pDmrResult_noHC_75 <- runFunctions(pDmr_noHC_75, p = 0.8, pos = "idiopathic_autism") # less 1, 286 predictors
pDmrResult_noHC_79 <- runFunctions(pDmr_noHC_79, p = 0.8, pos = "idiopathic_autism") # less 1, 286 predictors
pDmrResult_noHC_60 <- runFunctions(pDmr_noHC_60, p = 0.8, pos = "idiopathic_autism") # 270 predictors

# generate results table
rf_rett_hc90 <- sumRes(rDmrResult_noHC_90, caption = "0.90 highly correlated variables removed")
rf_rett_hc80 <- sumRes(rDmrResult_noHC_80, caption = "0.80 highly correlated variables removed")
rf_rett_hc70 <- sumRes(rDmrResult_noHC_70, caption = "0.70 highly correlated variables removed")
rf_rett_hc90
rf_rett_hc80
rf_rett_hc70

rf_dup_hc90 <- sumRes(dDmrResult_noHC_90, caption = "0.90 highly correlated variables removed")
rf_dup_hc80 <- sumRes(dDmrResult_noHC_80, caption = "0.80 highly correlated variables removed")
rf_dup_hc70 <- sumRes(dDmrResult_noHC_70, caption = "0.70 highly correlated variables removed")
rf_dup_hc90
rf_dup_hc80
rf_dup_hc70

rf_asd_hc90 <- sumRes(aDmrResult_noHC_90, caption = "0.90 highly correlated variables removed")
rf_asd_hc80 <- sumRes(aDmrResult_noHC_80, caption = "0.80 highly correlated variables removed")
rf_asd_hc70 <- sumRes(aDmrResult_noHC_70, caption = "0.70 highly correlated variables removed")
rf_asd_hc90
rf_asd_hc80
rf_asd_hc70 

rf_plac_hc70 <- sumRes(pDmrResult_noHC_70, caption = "0.70 to 0.79 highly correlated variables removed")
rf_plac_hc60 <- sumRes(pDmrResult_noHC_60, caption = "0.60 highly correlated variables removed")
rf_plac_hc70
rf_plac_hc60

# Neural Network Model Results --------------------------------------------
NNsumRes <- function(dmrResult, caption) {
  a <- paste("Model", dmrResult$nnModel$modelInfo$label, sep = ": ")
  b <- paste("Disorder", dmrResult$confMat$positive, sep = ": ")
  c <- paste(a, b, sep = ", ")
  ntrain <- dim(dmrResult$nnModel$trainingData)[1]
  ntest <- nrow(dmrResult$probPreds)
  npred <- dim(dmrResult$nnModel$trainingData)[2] - 1
  
  sumTable <- tibble(Measure = as.character(), Value = as.numeric()) %>%
    add_row(Measure = "Number of Samples in Training Data", Value = round(ntrain)) %>%
    add_row(Measure = "Number of Samples in Testing Data", Value = ntest) %>%
    add_row(Measure = "Number of Predictors", Value = npred) %>%
    add_row(Measure = "Accuracy", Value = dmrResult$confMat$overall["Accuracy"]) %>%
    add_row(Measure = "Kappa", Value = dmrResult$confMat$overall["Kappa"]) %>%
    add_row(Measure = "Accuracy P Value (Acc > NIR)", Value = dmrResult$confMat$overall["AccuracyPValue"]) %>%
    add_row(Measure = "Sensitivity", Value = dmrResult$confMat$byClass["Sensitivity"]) %>%
    add_row(Measure = "Specificity", Value = dmrResult$confMat$byClass["Specificity"]) %>%
    add_row(Measure = "Positive Predictive Values", Value = dmrResult$confMat$byClass["Pos Pred Value"]) %>%
    add_row(Measure = "Negative Predictive Values", Value = dmrResult$confMat$byClass["Neg Pred Value"]) %>%
    kable(caption = paste(c, caption, sep = " - ")) %>%
    kable_styling() %>%
    column_spec(1:2, color = "black") %>%
    add_header_above(header = c("Summarized results from classification algorithm" = 2), 
                     align = "c")
  return(sumTable)
}

NNrunFunctions(rDmr_noHC_90, p = 0.8, pos = "Rett") # didn't converge, too many weights warnings
nn_rDmrResult_noHC_80 <- NNrunFunctions(rDmr_noHC_80, p = 0.8, pos = "Rett") # converged, too many weights warnings
nn_rDmrResult_noHC_70 <- NNrunFunctions(rDmr_noHC_70, p = 0.8, pos = "Rett")

nn_dDmrResult_noHC_90 <- NNrunFunctions(dDmr_noHC_90, p = 0.8, pos = "Rett") # didn't converge, too many weights warnings
nn_dDmrResult_noHC_80 <- NNrunFunctions(dDmr_noHC_80, p = 0.8, pos = "Rett") # converged, error, subscript out of bounds
nn_dDmrResult_noHC_70 <- NNrunFunctions(dDmr_noHC_70, p = 0.8, pos = "Rett") # converged, error, subscript out of bounds

nn_aDmrResult_noHC_90 <- NNrunFunctions(aDmr_noHC_90, p = 0.8, pos = "Rett") # converged, error, subscript, weights warnings
nn_aDmrResult_noHC_80 <- NNrunFunctions(aDmr_noHC_80, p = 0.8, pos = "Rett") # converged, error, subscript, weights warnings
nn_aDmrResult_noHC_70 <- NNrunFunctions(aDmr_noHC_70, p = 0.8, pos = "Rett") # converged, error, subscript, weights warnings

nn_pDmrResult_noHC_70 <- NNrunFunctions(pDmr_noHC_70, p = 0.8, pos = "Rett") # converged, error, subscript, weights warnings
nn_pDmrResult_noHC_60 <- NNrunFunctions(pDmr_noHC_60, p = 0.8, pos = "Rett") # converged, error, subscript, weights warnings

nn_rett_hc80 <- NNsumRes(nn_rDmrResult_noHC_80, caption = "0.80 highly correlated variables removed, too many weights warnings")
nn_rett_hc70 <- NNsumRes(nn_rDmrResult_noHC_70, caption = "0.70 highly correlated variables removed")
nn_rett_hc80
nn_rett_hc70



# ran for more than 10 min
# Error: Stopping
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# > warnings()
# Warning messages:
#   1: model fit failed for Fold1.Rep01: size=1, decay=0e+00 Error in nnet.default(x, y, w, entropy = TRUE, ...) : 
#   too many (4644) weights
NNrunFunctions(rDmr, p = 0.8, pos = "Rett") # warnings, too many weights, data set too large

NNrunFunctions(rDmr_vi$sixty) # 10 samples, 28 predictors, size = 1, decay = 0.1, acc = 1, 1
NNrunFunctions(rDmr_vi$seventy) # 10 samples, 7 predictors, size = 1, decay = 0.1, acc = 1, 1
NNrunFunctions(rDmr_vi$eighty) # 10 samples, 7 predictors, size = 1, decay = 0.1, acc = 1, 1
NNrunFunctions(rDmr_vi$ninety) # 10 samples, 4 predictors, size = 1, decay = 0.1, acc = 1, 1

NNrunFunctions(dDmr_vi$sixty) # 8 samples, 32 predictors, size = 1, decay = 0.1, acc = 1, 1
NNrunFunctions(dDmr_vi$seventy) # 8 samples, 9 predictors, size = 1, decay = 0.1, acc = 1, 0.5
NNrunFunctions(dDmr_vi$eighty) # 8 samples, 9 predictors, size = 1, decay = 0.1, acc = 1, 0.5
NNrunFunctions(dDmr_vi$ninety) # 8 samples, 6 predictors, size = 1, decay = 0.1, acc = 1, 0.5

NNrunFunctions(aDmr_vi$sixty) # 22 samples, 4 predictors, size = 3, decay = 0.1, acc = 0.9083, 0.6
NNrunFunctions(aDmr_vi$seventy) # 22 samples, 1 predictors, size = 5, decay = 0.85059, acc = 1, 1
NNrunFunctions(aDmr_vi$eighty) # 22 samples, 1 predictors, size = 5, decay = 0.8505952, acc = 1, 1
NNrunFunctions(aDmr_vi$ninety) # 22 samples, 1 predictors, size = 1, decay = 0.8505952, acc = 1, 1


# Feature Selection - Variable Importance ---------------------------------

# run after selecting important variables
cutoff_vi <- c(60, 70, 80, 90)
rDmr_vi <- list() 
rDmr_vi$sixty <- selectImpVar(rDmr, rDmrResult$rfModel, cutoffValue = 60)
rDmr_vi$seventy <- selectImpVar(rDmr, rDmrResult$rfModel, cutoffValue = 70) # accuracy: 1 -> 0.5
rDmr_vi$eighty <- selectImpVar(rDmr, rDmrResult$rfModel, cutoffValue = 80) # accuracy: 1 -> 0.5
rDmr_vi$ninety <- selectImpVar(rDmr, rDmrResult$rfModel, cutoffValue = 90)

dDmr_vi <- list()
dDmr_vi$sixty <- selectImpVar(dDmr, dDmrResult$rfModel, cutoffValue = 60)
dDmr_vi$seventy <- selectImpVar(dDmr,dDmrResult$rfModel, cutoffValue = 70) # accuracy: 1 -> 1
dDmr_vi$eighty <- selectImpVar(dDmr, dDmrResult$rfModel, cutoffValue = 80) # accuracy: 1 -> 1
dDmr_vi$ninety <- selectImpVar(dDmr, dDmrResult$rfModel, cutoffValue = 90)

aDmr_vi <- list()
aDmr_vi$sixty <- selectImpVar(aDmr, aDmrResult$rfModel, cutoffValue = 60)
aDmr_vi$seventy <- selectImpVar(aDmr,aDmrResult$rfModel, cutoffValue = 70) # accuracy: 0.8 -> 1
aDmr_vi$eighty <- selectImpVar(aDmr, aDmrResult$rfModel, cutoffValue = 80) # accuracy: 0.8 -> 1
aDmr_vi$ninety <- selectImpVar(aDmr, aDmrResult$rfModel, cutoffValue = 90)

rf_rDmrResult_vi_60 <- runFunctions(rDmr_vi$sixty, p = 0.8, pos = "Rett")
rf_rDmrResult_vi_70 <- runFunctions(rDmr_vi$seventy, p = 0.8, pos = "Rett")
rf_rDmrResult_vi_80 <- runFunctions(rDmr_vi$eighty, p = 0.8, pos = "Rett")
rf_rDmrResult_vi_90 <- runFunctions(rDmr_vi$ninety, p = 0.8, pos = "Rett")

rf_dDmrResult_vi_60 <- runFunctions(dDmr_vi$sixty, p = 0.8, pos = "Dup15q")
rf_dDmrResult_vi_70 <- runFunctions(dDmr_vi$seventy, p = 0.8, pos = "Dup15q")
rf_dDmrResult_vi_80 <- runFunctions(dDmr_vi$eighty, p = 0.8, pos = "Dup15q")
rf_dDmrResult_vi_90 <- runFunctions(dDmr_vi$ninety, p = 0.8, pos = "Dup15q")

rf_aDmrResult_vi_60 <- runFunctions(aDmr_vi$sixty, p = 0.8, pos = "ASD")
rf_aDmrResult_vi_70 <- runFunctions(aDmr_vi$seventy, p = 0.8, pos = "ASD")
rf_aDmrResult_vi_80 <- runFunctions(aDmr_vi$eighty, p = 0.8, pos = "ASD")
rf_aDmrResult_vi_90 <- runFunctions(aDmr_vi$ninety, p = 0.8, pos = "ASD")



vi60pred <- colnames(rDmr_vi$sixty[-1])
rf_rett_vi60 <- sumRes(rf_rDmrResult_vi_60, caption = "variable importance 60")
rf_rett_vi70 <- sumRes(rf_rDmrResult_vi_70, caption = "variable importance 70")
rf_rett_vi80 <- sumRes(rf_rDmrResult_vi_80, caption = "variable importance 80")
rf_rett_vi90 <- sumRes(rf_rDmrResult_vi_90, caption = "variable importance 90")
rf_rett_vi60
rf_rett_vi70
rf_rett_vi80
rf_rett_vi90

rf_dup_vi60 <- sumRes(rf_dDmrResult_vi_60, caption = "variable importance 60")
rf_dup_vi70 <- sumRes(rf_dDmrResult_vi_70, caption = "variable importance 70")
rf_dup_vi80 <- sumRes(rf_dDmrResult_vi_80, caption = "variable importance 80")
rf_dup_vi90 <- sumRes(rf_dDmrResult_vi_90, caption = "variable importance 90")
rf_dup_vi60
rf_dup_vi70
rf_dup_vi80
rf_dup_vi90



rDmr_vi_Result <- lapply(rDmr_vi, runFunctions)
dDmr_vi_Result <- lapply(dDmr_vi, runFunctions)
# ∨ warnings() 50: In randomForest.default(x, y, mtry = param$mtry, ...) :invalid mtry: reset to within valid range
aDmr_vi_Result <- lapply(aDmr_vi, runFunctions)  
# aDmrResult_vi <- runFunctions(aDmr_vi$sixty), no error
# aDmrResult_vi <- runFunctions(aDmr_vi$seventy), eighty, ninety -> only 1 predictor
# mtry warning: Tuning parameter 'mtry' was held constant at a value of 2


# Compare models and accuracy before and after variable importance
rDmrResult$rfModel$results$Accuracy
rDmrResult$confMat$overall["Accuracy"]

# 28 predictors, accuracy = 1
rDmr_vi_Result$sixty$rfModel$results$Accuracy
rDmr_vi_Result$sixty$confMat$overall["Accuracy"]

# 7 predictors, accuracy = 1
rDmr_vi_Result$seventy$rfModel$results$Accuracy
rDmr_vi_Result$seventy$confMat$overall["Accuracy"]

# 7 predictors, accuracy = 1
rDmr_vi_Result$eighty$rfModel$results$Accuracy
rDmr_vi_Result$eighty$confMat$overall["Accuracy"]

# 4 predictors, accuracy = 1
rDmr_vi_Result$ninety$rfModel$results$Accuracy
rDmr_vi_Result$ninety$confMat$overall["Accuracy"]

# ROC curve, variable importance 
