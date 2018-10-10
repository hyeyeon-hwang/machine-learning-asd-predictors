# caret
# https://www.analyticsvidhya.com/blog/2016/12/practical-guide-to-implement-machine-learning-with-caret-package-in-r-with-practice-problem/

library(tidyverse)
library(caret)
library(randomForest)

rettDmrFull <- read.delim("Rett_sig_individual_smoothed_DMR_methylation.txt")
rettInfo <- read.csv("Rett_sample_info.csv")

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
    spread(seqId, values) %>% # rows to cols
    # add diagnosis column to front 
    add_column( diagnosis = factor(c("control", "control", "control", 
                  "rett", "rett", "rett", 
                  "rett", "rett", "rett", 
                  "control", "control", "control")), 
                .after = 1) 

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
vi_plot <- plot( vi, main = "Random Forest - Variable Importance")
vi_plot
