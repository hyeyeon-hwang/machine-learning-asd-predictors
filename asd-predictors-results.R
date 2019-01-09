source("asd-predictors.R")
library(knitr)
library(kableExtra)

# Random Forest Model Results ---------------------------------------------
rDmrResult <- runFunctions(rDmr)
plot.confusion(rDmrResult$confMat$table)

dDmrResult <- runFunctions(dDmr)
aDmrResult <- runFunctions(aDmr)
pDmrResult <- runFunctions(pDmr)
pDmrCBResult <- runFunctions(pDmrCB)

# tables using kabble
rDmrResult$confMat$table %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", 
                                      "hover",
                                      "condensed", 
                                      "responsive"),
                font_size = 15)

sumRes <- function(dmrResult, alg, disorder) {
  sumTable <- tibble(Measure = as.character(), Value = as.numeric()) %>%
    add_row(Measure = "Accuracy", Value = dmrResult$confMat$overall["Accuracy"]) %>%
    add_row(Measure = "Kappa", Value = dmrResult$confMat$overall["Kappa"]) %>%
    add_row(Measure = "Accuracy P Value (Acc > NIR)", Value = dmrResult$confMat$overall["AccuracyPValue"]) %>%
    add_row(Measure = "Sensitivity", Value = dmrResult$confMat$byClass["Sensitivity"]) %>%
    add_row(Measure = "Specificity", Value = dmrResult$confMat$byClass["Specificity"]) %>%
    add_row(Measure = "Positive Predictive Values", Value = dmrResult$confMat$byClass["Pos Pred Value"]) %>%
    add_row(Measure = "Negative Predictive Values", Value = dmrResult$confMat$byClass["Neg Pred Value"]) %>%
    kable() %>%
    kable_styling() %>%
    column_spec(1:2, color = "black") %>%
    add_header_above(header = c("Summarized results from classification algorithm" = 2), 
                     align = "c")
  return(sumTable)
}
rf_rett_table <- sumRes(rDmrResult, alg = "random forest", disorder = "rett")
rf_dup_table <- sumRes(dDmrResult, alg = "random forest", disorder = "dup")
rf_asd_table <- sumRes(aDmrResult, alg = "random forest", disorder = "asd")


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


# Feature Selection - Remove highly correlated variables ------------------

# run after removing highly correlated variables
rDmr_noHC_90 <- removeHighCor(rDmr, cutoffValue = 0.90) #10 samples, 2884 predictors, acc = 0.925, 1
dDmr_noHC_90 <- removeHighCor(dDmr, cutoffValue = 0.90) #8 samples, 1196 predictors, acc = 0.769, 0.5
aDmr_noHC_90 <- removeHighCor(aDmr, cutoffValue = 0.90) #22 samples, 470 predictors, acc = 0.849, 1

rDmr_noHC_80 <- removeHighCor(rDmr, cutoffValue = 0.80) #10 samples, 460 predictors, acc = 0.2027, 0
dDmr_noHC_80 <- removeHighCor(dDmr, cutoffValue = 0.80) #8 samples, 116 predictors, acc = 0.0694, 0
aDmr_noHC_80 <- removeHighCor(aDmr, cutoffValue = 0.80) #22 samples, 443 predictors, acc = 0.87, 1

# 0.75 cutoff: all models drastically worsen
# 0.85 cutoff: only asd better 
# 0.90 cutoff: all models slightly better
rDmrResult_noHC_90 <- runFunctions(rDmr_noHC_90) 
dDmrResult_noHC_90 <- runFunctions(dDmr_noHC_90)
aDmrResult_noHC_90 <- runFunctions(aDmr_noHC_90)

rDmrResult_noHC_80 <- runFunctions(rDmr_noHC_80) 
dDmrResult_noHC_80 <- runFunctions(dDmr_noHC_80)
aDmrResult_noHC_80 <- runFunctions(aDmr_noHC_80)





# Neural Network Model Results --------------------------------------------

# ran for more than 10 min
# Error: Stopping
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# > warnings()
# Warning messages:
#   1: model fit failed for Fold1.Rep01: size=1, decay=0e+00 Error in nnet.default(x, y, w, entropy = TRUE, ...) : 
#   too many (4644) weights
NNrunFunctions(rDmr) # warnings, too many weights, data set too large

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



