library(tidyverse)
library(caret)
library(randomForest)

rettDmrFull <- read.delim("data/Rett_sig_individual_smoothed_DMR_methylation.txt", check.names = FALSE)
rettInfo <- read.csv("data/Rett_sample_info.csv")

dupDmrFull <- read.delim("data/Dup15q_sig_individual_smoothed_DMR_methylation.txt")
dupInfo <- read.csv("data/Dup15q_sample_info.csv")

asdDmrFull <- read.delim("data/ASD_sig_individual_smoothed_DMR_methylation.txt")
asdInfo <- read.csv("data/ASD_sample_info.csv")

# DMR Dataset -------------------------------------------------------------
# exclude range of columns from 'width' to 'RawDiff', transpose, add diagnosis column
dmr <- rettDmrFull %>% 
    as.tibble() %>% 
    select(-(width:RawDiff)) %>%
    unite(seqId1, seqnames, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") %>%
    # transpose: cols to rows
    gather(sampleID, values, -seqId) %>% # cols to rows
    # transpose: rows to cols
    spread(seqId, values) 

# list is like a dictionary
# add diagnosis column, match and read in from file
dmr %>%
    add_column( diagnosis = "-", .after = 1 )

# %>% # rows to cols
#     # add diagnosis column to front, read in from file later
#     add_column( diagnosis = factor(c("control", "control", "control",
#                   "rett", "rett", "rett",
#                   "rett", "rett", "rett",
#                   "control", "control", "control")),
#                 .after = 1)

# remove sampleID column: 1136, 1406, 1711, 1815, 4687, 4852, 
#                         5020, 5075, 5214, 738, 754, 812
rettDmr <- dmr %>%
    select(-sampleID)

# small DMR data set 
rettDmr <- rettDmr[, 1:100]

# Partition data into training and testing --------------------------------
seed <- 9999
set.seed(seed)
trainIndex <- createDataPartition(rettDmr$diagnosis, 
                                  p = 0.8,
                                  list = FALSE )

training <- rettDmr[trainIndex, ]
testing <- rettDmr[-trainIndex, ]

# Model: Random Forest different trControl  ---------------------------
fitControl <- trainControl(method = "none", returnResamp = "final")
set.seed(seed)
rf_default <- train( diagnosis ~ ., 
                     data = training, 
                     method = "rf", 
                     trControl = fitControl)

# Model: Random Forest 10-fold Cross Validation ---------------------------
fitControl <- trainControl(method = "repeatedcv", # 10-fold cv
                           number = 10, 
                           repeats = 10, 
                           classProbs = TRUE) 
set.seed(seed)
rf_default <- train( diagnosis ~ ., 
                     data = training, 
                     method = "rf", 
                     trControl = fitControl)

# predict the outcome on a test set
rettPredict <- predict(rf_default, testing)
confusionMatrix(rettPredict, testing$diagnosis)

# Model: Random Forest using ranger ---------------------------------------
library(ranger)
rf_fit <- train(diagnosis ~ ., 
                data = training, 
                method = "ranger")
rf_fit

# predict the outcome on a test set
rangerPredict <- predict(rf_fit, testing)
confusionMatrix(rangerPredict, testing$diagnosis)

# Model: Stochastic Gradient Boosting -------------------------------------
library(gbm)
gbmFitControl <- trainControl(method = "repeatedcv", 
                              number = 2, 
                              repeats = 2)
set.seed(seed)
gbmFit1 <- train(diagnosis ~ ., 
                 data = training, 
                 method = "gbm", 
                 trControl = gbmFitControl, 
                 verbose = FALSE)

# Variable Importance for Random Forest -----------------------------------
vi <- varImp(object = rf_default)
vi
vi_plot <- plot(vi, main = "Random Forest - Variable Importance")
vi_plot

# warnings test: no warnings if no resampling -----------------------------------------------------------
# below fitControl and train give no resampling, no warnings
# but why no resampling? 
# warning ok? https://github.com/topepo/caret/issues/905

fitControl <- trainControl(method = "none", 
                           classProbs = TRUE, 
                           returnData = TRUE, 
                           returnResamp = "all", 
                           savePredictions = "all") 
fitControl

set.seed(seed)
rf_noresamp <- train( diagnosis ~ .,
                      data = training, 
                      method = "rf", 
                      trControl = fitControl) 
rf_noresamp 

noresampPredict <- predict(rf_noresamp, testing) 
confusionMatrix(noresampPredict, testing$diagnosis) 

# warnings: 
# 1. missing values in resampled performance measures
# caret/workflows.R
# resampleHist.R
# postResample.R
# resampleSummary.R
# resamples. R - collation and visualization of resampling results 

# caret:
# https://www.analyticsvidhya.com/blog/2016/12/practical-guide-to-implement-machine-learning-with-caret-package-in-r-with-practice-problem/