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

rfModel <- list()
rfModel$rett <- fitRfModel(dmr$rett)
rfModel$dup <- fitRfModel(dmr$dup)
rfModel$asd <- fitRfModel(dmr$asd)
rfModel$plac <- fitRfModel(dmr$plac)

rfModel$rett
rfModel$rett$finalModel
rfModel$rett$pred # savePredictions = "final" outputs predicted probabilites for resamples with optimal mtry
confusionMatrix.train(rfModel$rett)

rfModel$dup
rfModel$dup$finalModel
rfModel$dup$pred
confusionMatrix.train(rfModel$dup)

rfModel$asd 
rfModel$asd$finalModel
rfModel$asd$pred
confusionMatrix.train(rfModel$asd)

rfModel$plac
rfModel$plac$finalModel
rfModel$plac$pred
confusionMatrix.train(rfModel$plac)


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