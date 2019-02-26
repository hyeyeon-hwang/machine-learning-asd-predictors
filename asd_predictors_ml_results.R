source("asd_predictors_ml.R") #start 11:38 - 11:44
library(knitr)
library(kableExtra)


# cmTables ----------------------------------------------------------------

cmTable <- function(df, colNum) {
  df %>%
    kable() %>%
    kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
    add_header_above(c(" ", "Reference" = colNum)) %>%
    add_header_above(header = c("5-fold Cross Validated Confusion Matrix" = colNum + 1), align = "c") %>%
    row_spec(0, bold = F) # still bolded
}
cmPlacFile <- cmTable(cmPlacDf, colNum = 2) 
cmAsdFile <- cmTable(cmAsdCnt$table, colNum = 2)
cmDupFile <- cmTable(cmDupDf, colNum = 2)
cmRettFile <- cmTable(cmRettDf, colNum = 2)
cmMi3File <- cmTable(cmMi3Df, colNum = 4)
cmMi3gFile <- cmTable(cmMi3gDf, colNum = 2)
cmMi4File <- cmTable(cmMi4Df, colNum = 5)
cmMi4gFile <- cmTable(cmMi4gDf, colNum = 2)

save_kable(cmPlacFile, "../result_figures/cmPlacFile.pdf")
save_kable(cmAsdFile, "../result_figures/cmAsdFile.pdf")
save_kable(cmDupFile, "../result_figures/cmDupFile.pdf")
save_kable(cmRettFile, "../result_figures/cmRettFile.pdf")
save_kable(cmMi3File, "../result_figures/cmMi3File.pdf")
save_kable(cmMi3gFile, "../result_figures/cmMi3gFile.pdf")
save_kable(cmMi4File, "../result_figures/cmMi4File.pdf")
save_kable(cmMi4gFile, "../result_figures/cmMi4gFile.pdf")

# add accuracy and kappa values to this table for each fold?
# probRett ----------------------------------------------------------------
probRett <- rfModel$rett$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, Rett, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_Rett = Rett, 
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  #collapse_rows(columns = 1, valign = "top")
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Rett Dataset" = 6), align = "c")

# probDup -----------------------------------------------------------------
probDup <- rfModel$dup$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, Dup15q, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_Dup15q = Dup15q, 
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Dup15q Dataset " = 6), align = "c")

# probAsd -----------------------------------------------------------------
probAsd <- rfModel$asd$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, ASD, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_ASD = ASD, 
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for ASD Dataset " = 6), align = "c")


# probPlac ----------------------------------------------------------------
probPlac <- rfModel$plac$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, Idiopathic_autism, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_Idiopathic_autism = Idiopathic_autism, 
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Placenta Dataset " = 6), align = "c")


# probMi3 -----------------------------------------------------------------
probMi3 <- rfModel$mi3$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, ASD, Dup15q, Rett, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_ASD = ASD, 
         Probability_Dup15q = Dup15q,
         Probability_Rett = Rett,
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Merged ASD-Dup15q-Rett Dataset" = 8), align = "c")

# probMi3g -----------------------------------------------------------------
probMi3g <- rfModel$mi3Grouped$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, Positive, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_Positive = Positive, 
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Merged ASD-Dup15q-Rett Dataset (with 2 classes)" = 6), align = "c")

# probMi4 -----------------------------------------------------------------
probMi4 <- rfModel$mi4$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, ASD, Dup15q, Rett, Idiopathic_autism, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_ASD = ASD, 
         Probability_Dup15q = Dup15q,
         Probability_Rett = Rett,
         Probability_Idiopathic_autism = Idiopathic_autism,
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Merged ASD-Dup15q-Rett-Placenta Dataset" = 9), align = "c") 

# probMi4g ----------------------------------------------------------------
probMi4g <- rfModel$mi4Grouped$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, Positive, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_Positive = Positive, 
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Merged ASD-Dup15q-Rett-Placenta Dataset (with 2 classes)" = 6), align = "c")









# save prob to pdf -------------------------------------------------------------
save_kable(probRett, "../result_figures/probRett.pdf")
save_kable(probDup, "../result_figures/probDup.pdf")
save_kable(probAsd, "../result_figures/probAsd.pdf")
save_kable(probPlac, "../result_figures/probPlac.pdf")
save_kable(probMi3, "../result_figures/probMi3.pdf")
save_kable(probMi3g, "../result_figures/probMi3g.pdf")
save_kable(probMi4, "../result_figures/probMi4.pdf")
save_kable(probMi4g, "../result_figures/probMi4g.pdf")


# sumMi4 ------------------------------------------------------------------

nsamp <- dim(rfModel$mi4$trainingData)[1]
npred <- dim(rfModel$mi4$trainingData)[2] - 1
tn <- cmMi4Cnt$table["Control", "Control"]
tp <- cmMi4Cnt$table["ASD", "ASD"] + 
  cmMi4Cnt$table["Dup15q", "Dup15q"] +
  cmMi4Cnt$table["Rett", "Rett"] + 
  cmMi4Cnt$table["Idiopathic_autism", "Idiopathic_autism"]
fn <- cmMi4Cnt$table["Control", "ASD"] + 
  cmMi4Cnt$table["Control", "Dup15q"] +
  cmMi4Cnt$table["Control", "Rett"] +
  cmMi4Cnt$table["Control", "Idiopathic_autism"] 
fp <- cmMi4Cnt$table["ASD", "Control"] +
  cmMi4Cnt$table["Dup15q", "Control"] +
  cmMi4Cnt$table["Rett", "Control"] +
  cmMi4Cnt$table["Idiopathic_autism", "Control"]
acc <- (tn + tp) / (tn + tp + fn + fp)
sens <- tp / (tp + fn)
spec <- tn / (tn + fp)
ppv <- tp / (tp + fp)
npv <- tn / (tn + fn)

# sumPlac ------------------------------------------------
nsamp <- dim(rfModel$plac$trainingData)[1]
npred <- dim(rfModel$plac$trainingData)[2] - 1
tn <- cmPlacCnt$table["Control","Control"] 
tp <- cmPlacCnt$table["Idiopathic_autism", "Idiopathic_autism"]
fn <- cmPlacCnt$table["Control","Idiopathic_autism"]
fp <- cmPlacCnt$table["Idiopathic_autism", "Control"]
acc <- (tn + tp) / (tn + tp + fn + fp)
sens <- tp / (tp + fn)
spec <- tn / (tn + fp)
ppv <- tp / (tp + fp)
npv <- tn / (tn + fn)

sumPlac <- tibble(Measure = as.character(), Value = as.numeric()) %>%
  add_row(Measure = "Total number of samples", Value = nsamp) %>%
  add_row(Measure = "Total number of predictors", Value = npred) %>%
  add_row(Measure = "Accuracy", Value = acc) %>%
  add_row(Measure = "Sensitivity", Value = sens) %>%
  add_row(Measure = "Specificity", Value = spec) %>%
  add_row(Measure = "Positive predictive value", Value = ppv) %>%
  add_row(Measure = "Negative predictive value", Value = npv) %>%
  kable(caption = paste(c, caption, sep = " - ")) %>%
  kable_styling() %>%
  column_spec(1:2, color = "black") %>%
  add_header_above(header = c("Summarized results from classification algorithm" = 2), 
                   align = "c")
  

sumRes <- function(dmrResult, colNum, caption) {
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
    column_spec(1:colNum, color = "black") %>%
    add_header_above(header = c("Summarized results from classification algorithm" = colNum+1), 
                     align = "c")
  return(sumTable)
}


# in progress -------------------------------------------------------------
save_kable(probMi4g, "probMi4g.pdf")

  

