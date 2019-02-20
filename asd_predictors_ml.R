source("asd_predictors_data_prep.R") # ~1min

seed <- 9999
# 5 fold cross validation -------------------------------------------------
library(caret)
fitControl_5fold <- trainControl(method = "cv", 
                                 number = 5, 
                                 verboseIter = TRUE,
                                 returnResamp = "final",#all
                                 savePredictions = "final", #"all"
                                 classProbs = TRUE) 

fitRfModel <- function(dmrData) {
  set.seed(seed)
  model <- train(diagnosis ~ .,
                 data = dmrData[,-1],
                 method = "rf",
                 trControl = fitControl_5fold)
  return(model)
}

fitRfModelVita <- function(dmrData) {
  set.seed(seed)
  model <- train(diagnosis ~ .,
                 data = dmrData[,-1],
                 method = "rf",
                 trControl = fitControl_5fold)
  return(model)
}
# https://stackoverflow.com/questions/33470373/applying-k-fold-cross-validation-model-using-caret-package
# when you perform k-fold cv, you are already making a
# prediction for each sample, just over 5 diff models (k = 5)
# so there is no need to make a prediction on the complete data,
# as you already have their predictions from the k diff models


# https://stats.stackexchange.com/questions/52274/how-to-choose-a-predictive-model-after-k-fold-cross-validation
# the purpose of cv is not to come up with final model, we don't do any real prediction
# for that, we want to use all the data to build the best model
# purpose of cv: model checking, not model building 
# use cv to choose better performing model 
#     then train new model on all the data -> then predict on new data
  

rfModel <- list()
# new rett input has X's before sample name, sampleName doesn't match
# so NA for diagnosis
# changed new rett input file to exclude the X's
rfModel$rett <- fitRfModel(dmr$rett) # 2  1         1  
rfModel$dup <- fitRfModel(dmr$dup) # 2  1         1 
rfModel$asd <- fitRfModel(dmr$asd)  # 2   0.9266667  0.8233766
rfModel$plac <- fitRfModel(dmr$plac) # 2   0.9750000  0.9500000
rfModel$mi3 <- fitRfModel(dmr$mi3)
rfModel$mi3Grouped <- fitRfModel(dmr$mi3Grouped)
rfModel$mi4 <- fitRfModel(dmr$mi4)
rfModel$mi4Grouped <- fitRfModel(dmr$mi4Grouped)

# finalModel calls randomForest() with optimal mtry
# my results: from confusionMatrix.train() = results from resamples, in percentages


rfModel$mi4Grouped # mtry, acc, kappa: 2 0.8156863  0.2260369
rfModel$mi4Grouped$finalModel
rfModel$mi4Grouped$pred
cmMi4gPerc <- confusionMatrix.train(rfModel$mi4Grouped) 
cmMi4gCnt <- confusionMatrix.train(rfModel$mi4Grouped, norm = "none")

rfModel$mi4 # mtry, acc, kappa: 2 0.4906046  0.3584639
rfModel$mi4$finalModel
rfModel$mi4$pred
cmMi4Perc <- confusionMatrix.train(rfModel$mi4) 
cmMi4Cnt <- confusionMatrix.train(rfModel$mi4, norm = "none")

rfModel$mi3Grouped #115  0.7233333  0.3500248
rfModel$mi3Grouped$finalModel
rfModel$mi3Grouped$pred
cmMi3gPerc <- confusionMatrix.train(rfModel$mi3Grouped) # don't need?
cmMi3gCnt <- confusionMatrix.train(rfModel$mi3Grouped, norm = "none")

rfModel$mi3 # 115  0.6250000  0.4137285
rfModel$mi3$finalModel
rfModel$mi3$pred
cmMi3Perc <- confusionMatrix.train(rfModel$mi3)
cmMi3Cnt <- confusionMatrix.train(rfModel$mi3, norm = "none")

rfModel$rett # 2  1         1
rfModel$rett$finalModel
rfModel$rett$pred # savePredictions = "final" outputs predicted probabilites for resamples with optimal mtry
cmRettPerc <- confusionMatrix.train(rfModel$rett)
cmRettCnt <- confusionMatrix.train(rfModel$rett, norm = "none")

rfModel$dup # 2  1         1 
rfModel$dup$finalModel
rfModel$dup$pred
cmDupPerc <- confusionMatrix.train(rfModel$dup)
cmDupCnt <- confusionMatrix.train(rfModel$dup, norm = "none")

rfModel$asd # 2   0.9266667  0.8233766
rfModel$asd$finalModel
rfModel$asd$pred
cmAsdPerc <- confusionMatrix.train(rfModel$asd)
cmAsdCnt <- confusionMatrix.train(rfModel$asd, norm = "none")

rfModel$plac # 2   0.9750000  0.9500000
rfModel$plac$finalModel
rfModel$plac$pred
cmPlacPerc <- confusionMatrix.train(rfModel$plac)
cmPlacCnt <- confusionMatrix.train(rfModel$plac, norm = "none", positive = "idiopathic_autism")


# placenta
cmPlacDf <- data.frame(Idiopathic_autism = c(cmPlacCnt$table[,"idiopathic_autism"][[1]], cmPlacCnt$table[,"idiopathic_autism"][[2]]),
                 Control = c(cmPlacCnt$table[,"control"][[1]], cmPlacCnt$table[,"control"][[2]]))
row.names(cmPlacDf) <- c("Idiopathic_autism", "Control")
# dup
cmDupDf <- data.frame(Dup15q = c(cmDupCnt$table[,"Dup15q"][[1]], cmDupCnt$table[,"Dup15q"][[2]]),
                       Control = c(cmDupCnt$table[,"Control"][[1]], cmDupCnt$table[,"Control"][[2]]))
row.names(cmDupDf) <- c("Dup15q", "Control")
# rett
cmRettDf <- data.frame(Rett = c(cmRettCnt$table[,"Rett"][[1]], cmRettCnt$table[,"Rett"][[2]]),
                      Control = c(cmRettCnt$table[,"Control"][[1]], cmRettCnt$table[,"Control"][[2]]))
row.names(cmRettDf) <- c("Rett", "Control")
# mi3Mi3
cmMi3Df <- data.frame(ASD = c(cmMi3Cnt$table[, "ASD"][[1]],
                              cmMi3Cnt$table[, "ASD"][[2]], 
                              cmMi3Cnt$table[, "ASD"][[3]],
                              cmMi3Cnt$table[, "ASD"][[4]]),
                      Dup15q = c(cmMi3Cnt$table[, "Dup15q"][[1]],
                                 cmMi3Cnt$table[, "Dup15q"][[2]],
                                 cmMi3Cnt$table[, "Dup15q"][[3]],
                                 cmMi3Cnt$table[, "Dup15q"][[4]]),
                      Rett = c(cmMi3Cnt$table[, "Rett"][[1]],
                               cmMi3Cnt$table[, "Rett"][[2]],
                               cmMi3Cnt$table[, "Rett"][[3]],
                               cmMi3Cnt$table[, "Rett"][[4]]),
                      Control = c(cmMi3Cnt$table[, "Control"][[1]],
                                  cmMi3Cnt$table[, "Control"][[2]],
                                  cmMi3Cnt$table[, "Control"][[3]],
                                  cmMi3Cnt$table[, "Control"][[4]]))
row.names(cmMi3Df) <- c("ASD", "Dup15q", "Rett", "Control")

library(kableExtra)
cmTable <- function(df, colNum) {
  # if(col == 4) {
  #   print("a")
  #   cmMi3Df %>%
  #     kable() %>%
  #     kable_styling(font_size = 12, full_width = F) %>%
  #     add_header_above(c(" ", "Reference" = colNum)) %>%
  #     add_header_above(header = c("5-fold Cross Validated Confusion Matrix" = colNum + 1), align = "c") 
  # }
  df %>%
    kable() %>%
    kable_styling(font_size = 12, full_width = F) %>%
    add_header_above(c(" ", "Reference" = colNum)) %>%
    add_header_above(header = c("5-fold Cross Validated Confusion Matrix" = colNum + 1), align = "c") 
}
cmTable(cmPlacDf, colNum = 2) 
cmTable(cmAsdCnt$table, colNum = 2)
cmTable(cmDupDf, colNum = 2)
cmTable(cmRettDf, colNum = 2)
cmTable(cmMi3Df, col = 4)



cmPlacCnt$table %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12, full_width = F) %>%
  #add_header_above(c(" ", "Reference" = 2)) %>%
  collapse_rows(columns = 1:2) %>%
  add_header_above(header = c("5-fold Cross Validated Confusion Matrix" = 3), align = "c") 
# Feature Selection - boruta ----------------------------------------------
library(Boruta)
set.seed(seed)
boruta.bank_train <- Boruta()


# Feature Selection - Vita ------------------------------------------------
library(vita)
#vitaRettRf <- rfModel$rett
#class(vitaRettRf) <- "randomForest"
vitaRett <- compVarImp(X = xRett, y = yRett, rForest = vitaRettRf, nPerm = 0)

#x = minus first 2 col
#y = 2nd column
xRett <- dmr$rett[, -(1:2)]
yRett <- dmr$rett[, 2]



# dmr$cb testing ----------------------------------------------------------
# dim(dmr$cb) = 46 97610, dim(dmr$rett) = 12 4643
rfModel$cb <- fitRfModel(dmr$cb) #start 3:17
# error: vector memory exhausted (limit reached?)
# 10 sampled rows -> error: vector memory exhausted 
dmrCbHc70 <- removeHighCor(dmr$cb, 70) 
# error: vector mem exh when making corMatrix :(

# cb 10k predictors / 97,608
set.seed(seed)
colIdx <- sample(1:ncol(dmr$cb), 10000)
cb10k <- dmr$cb[, colIdx] %>% 
  as.tibble() %>%
  add_column(sampleID = dmr$cb$sampleID, .before = 1) %>%
  add_column(diagnosis = dmr$cb$diagnosis, .after = 1)

rfModel$cb10k <- fitRfModel(cb10k) #start 4:00 - < 4:02
rfModel$cb10k$finalModel
rfModel$cb10k$pred
confusionMatrix.train(rfModel$cb10k)

# cb 20k predictors / 97,608 => ErrorL protect(): protecton stack overflow
# cb 15k predictors / 97,608
set.seed(seed)
colIdx2 <- sample(1:ncol(dmr$cb), 15000)
cb15k <- dmr$cb[, colIdx2] %>% 
  as.tibble() %>%
  add_column(sampleID = dmr$cb$sampleID, .before = 1) %>%
  add_column(diagnosis = dmr$cb$diagnosis, .after = 1)

rfModel$cb15k <- fitRfModel(cb15k) #start 4:07 - 4:10
rfModel$cb15k
rfModel$cb15k$finalModel
rfModel$cb15k$pred
confusionMatrix.train(rfModel$cb15k)


#  dmr$cb with classes: Positive, Control (not ASD, Control, Dup,  --------
dmr$cb2 <- dmr$cb 


# predictors/DMRs that overlap btwn brain & placenta samples --------------


  
# below old ---------------------------------------------------------------

predModel <- predict(rfModel$asd, aDmr)
confusionMatrix(predModel, aDmr$diagnosis) #complete overfitting




# obtained final random forest model for rett
# combined confusion matrix for resamples
# predict on rett_dDmr (make fake dataset that labels "dup15q" as "rett")
rett_dDmr <- dDmr %>% 
  mutate(diagnosis = str_replace(diagnosis, "Dup15q", "Rett"))
predict(rett_rfModel, rett_dDmr)

finalModel_asd <- asd_rfModel$finalModel$votes[,1]
finalModel_cont <- asd_rfModel$finalModel$votes[,2]
order <- asd_rfModel$pred[order(asd_rfModel$pred$rowIndex),]
pred_asd <- order$ASD
pred_cont <- order$Control

df_asd <- data.frame("finalModel_asd" = finalModel_asd, "pred_asd" = pred_asd)
