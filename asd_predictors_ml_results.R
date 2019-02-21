source("asd_predictors_ml.R") #start 11:38 - 11:44
library(knitr)
library(kableExtra)


# cmTables ----------------------------------------------------------------

cmTable <- function(df, colNum) {
  df %>%
    kable() %>%
    kable_styling(font_size = 12, full_width = F) %>%
    add_header_above(c(" ", "Reference" = colNum)) %>%
    add_header_above(header = c("5-fold Cross Validated Confusion Matrix" = colNum + 1), align = "c") %>%
    row_spec(0, bold = F) # still bolded
}
cmTable(cmPlacDf, colNum = 2) 
cmTable(cmAsdCnt$table, colNum = 2)
cmTable(cmDupDf, colNum = 2)
cmTable(cmRettDf, colNum = 2)
cmTable(cmMi3Df, colNum = 4)
cmTable(cmMi3gDf, colNum = 2)
cmTable(cmMi4Df, colNum = 5)
cmTable(cmMi4gDf, colNum = 2)

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
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
  #collapse_rows(columns = 1, valign = "top")
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Rett Dataset" = 6), align = "c")

# probDup -----------------------------------------------------------------
probDup <- rfModel$dup$pred %>%
  as.tibble() %>%
  arrange(Resample, rowIndex) %>%
  select(Resample, rowIndex, Dup15q, Control, pred, obs) %>%
  rename(Fold = Resample, 
         Sample_number = rowIndex, 
         Probability_Dup = Dup15q, 
         Probability_Control = Control, 
         Predicted_class = pred, 
         Actual_class = obs) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
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
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
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
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
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
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
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
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
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
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
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
  kable_styling(bootstrap_options = c("condensed"), font_size = 12) %>%
  add_header_above(header = c("Predicted Probabilities and Classes per Fold for Merged ASD-Dup15q-Rett-Placenta Dataset (with 2 classes)" = 6), align = "c")





# summary results for Plac ------------------------------------------------
nsamp <- dim(rfModel$plac$trainingData)[1]
npred <- dim(rfModel$Plac$trainingData)[2] - 1
tn <- cmPlacCnt$table["Control","Control"] 
tp <- cmPlacCnt$table["Idiopathic_autism", "Idiopathic_autism"]
fn <- cmPlacCnt$table["Control","Idiopathic_autism"]
fp <- cmPlacCnt$table["Idiopathic_autism", "Control"]
acc <- (tn + tp) / (tn + tp + fn + fp)

sumPlac <- tibble(Measure = as.character(), Value = as.numeric()) %>%
  add_row(Measure = "Total number of samples", Value = nsamp) %>%
  add_row(Measure = "Total number of predictors", Value = npred) %>%
  add_row(Measure = "Accuracy", Value = acc)


# in progress -------------------------------------------------------------
save_kable(probMi4g, "probMi4g.pdf")

  

