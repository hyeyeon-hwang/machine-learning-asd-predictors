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

# sumRes function ---------------------------------------------------------
sumRes <- function(dataType) {
  acc <- (tn + tp) / (tn + tp + fn + fp)
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  sumResTable <- tibble(Measure = as.character(), Value = as.numeric()) %>%
    add_row(Measure = "Total number of samples", Value = nsamp) %>%
    add_row(Measure = "Total number of predictors", Value = npred) %>%
    add_row(Measure = "Accuracy", Value = acc) %>%
    add_row(Measure = "Sensitivity", Value = sens) %>%
    add_row(Measure = "Specificity", Value = spec) %>%
    add_row(Measure = "Positive predictive value", Value = ppv) %>%
    add_row(Measure = "Negative predictive value", Value = npv) %>%
    kable() %>%
    kable_styling(bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
    column_spec(1:2, color = "black") %>%
    add_header_above(header = sumHeaderDesc(dataType), align = "c") 
  return(sumResTable)
}


sumHeaderDesc <- function(dataType) {
  if(dataType == "mi4") {
    header <- c("Summarized results from random forest algorithm on\n combined Asd-Dup15q-Rett-Placenta dataset" = 2)
  }
  if(dataType == "mi4g") {
    header <- c("Summarized results from random forest algorithm on\n combined Asd-Dup15q-Rett-Placenta dataset (with 2 classes)" = 2)
  }
  if(dataType == "mi3") {
    header <- c("Summarized results from random forest algorithm on\n combined Asd-Dup15q-Rett dataset" = 2)
  }
  if(dataType == "mi3g") {
    header <- c("Summarized results from random forest algorithm on\n combined Asd-Dup15q-Rett dataset (with 2 classes)" = 2)
  }
  if(dataType == "asd") {
    header <- c("Summarized results from random forest algorithm on\n ASD dataset" = 2)
  }
  if(dataType == "dup15q") {
    header <- c("Summarized results from random forest algorithm on\n Dup15q dataset" = 2)
  }
  if(dataType == "rett") {
  header <- c("Summarized results from random forest algorithm on\n Rett dataset" = 2)
  }
  if(dataType == "plac") {
  header <- c("Summarized results from random forest algorithm on\n Placenta dataset" = 2)
  }
  return(header)
}

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
sumMi4 <- sumRes("mi4")

# sumMi3 ------------------------------------------------------------------
nsamp <- dim(rfModel$mi3$trainingData)[1]
npred <- dim(rfModel$mi3$trainingData)[2] - 1
tn <- cmMi3Cnt$table["Control", "Control"]
tp <- cmMi3Cnt$table["ASD", "ASD"] + 
  cmMi3Cnt$table["Dup15q", "Dup15q"] +
  cmMi3Cnt$table["Rett", "Rett"]
fn <- cmMi3Cnt$table["Control", "ASD"] + 
  cmMi3Cnt$table["Control", "Dup15q"] +
  cmMi3Cnt$table["Control", "Rett"]
fp <- cmMi3Cnt$table["ASD", "Control"] +
  cmMi3Cnt$table["Dup15q", "Control"] +
  cmMi3Cnt$table["Rett", "Control"]
sumMi3 <- sumRes("mi3")

# sumMi4g ------------------------------------------------
nsamp <- dim(rfModel$mi4Grouped$trainingData)[1]
npred <- dim(rfModel$mi4Grouped$trainingData)[2] - 1
tn <- cmMi4gCnt$table["Control","Control"] 
tp <- cmMi4gCnt$table["Positive", "Positive"]
fn <- cmMi4gCnt$table["Control","Positive"]
fp <- cmMi4gCnt$table["Positive", "Control"]
sumMi4g <- sumRes("mi4g")

# sumMi3g ------------------------------------------------
nsamp <- dim(rfModel$mi3Grouped$trainingData)[1]
npred <- dim(rfModel$mi3Grouped$trainingData)[2] - 1
tn <- cmMi3gCnt$table["Control","Control"] 
tp <- cmMi3gCnt$table["Positive", "Positive"]
fn <- cmMi3gCnt$table["Control","Positive"]
fp <- cmMi3gCnt$table["Positive", "Control"]
sumMi3g <- sumRes("mi3g")

# sumAsd ------------------------------------------------
nsamp <- dim(rfModel$asd$trainingData)[1]
npred <- dim(rfModel$asd$trainingData)[2] - 1
tn <- cmAsdCnt$table["Control","Control"] 
tp <- cmAsdCnt$table["ASD", "ASD"]
fn <- cmAsdCnt$table["Control","ASD"]
fp <- cmAsdCnt$table["ASD", "Control"]
sumAsd <- sumRes("asd") 

# sumDup15q ------------------------------------------------
nsamp <- dim(rfModel$dup$trainingData)[1]
npred <- dim(rfModel$dup$trainingData)[2] - 1
tn <- cmDupCnt$table["Control","Control"] 
tp <- cmDupCnt$table["Dup15q", "Dup15q"]
fn <- cmDupCnt$table["Control","Dup15q"]
fp <- cmDupCnt$table["Dup15q", "Control"]
sumDup <- sumRes("dup15q") 

# sumRett ------------------------------------------------
nsamp <- dim(rfModel$rett$trainingData)[1]
npred <- dim(rfModel$rett$trainingData)[2] - 1
tn <- cmRettCnt$table["Control","Control"] 
tp <- cmRettCnt$table["Rett", "Rett"]
fn <- cmRettCnt$table["Control","Rett"]
fp <- cmRettCnt$table["Rett", "Control"]
sumRett <- sumRes("rett")

# sumPlac ------------------------------------------------
nsamp <- dim(rfModel$plac$trainingData)[1]
npred <- dim(rfModel$plac$trainingData)[2] - 1
tn <- cmPlacCnt$table["Control","Control"] 
tp <- cmPlacCnt$table["Idiopathic_autism", "Idiopathic_autism"]
fn <- cmPlacCnt$table["Control","Idiopathic_autism"]
fp <- cmPlacCnt$table["Idiopathic_autism", "Control"]
sumPlac <- sumRes("plac")

# save sumRes as pdf -------------------------------------------------------------
save_kable(sumRett, "../result_figures/sumRett.pdf")
save_kable(sumDup, "../result_figures/sumDup.pdf")
save_kable(sumAsd, "../result_figures/sumAsd.pdf")
save_kable(sumPlac, "../result_figures/sumPlac.pdf")
save_kable(sumMi3, "../result_figures/sumMi3.pdf")
save_kable(sumMi3g, "../result_figures/sumMi3g.pdf")
save_kable(sumMi4, "../result_figures/sumMi4.pdf") #, 
            #vwidth = "3cm", vheight = "3cm")
save_kable(sumMi4g, "../result_figures/sumMi4g.pdf")

  

