# caret
# https://www.analyticsvidhya.com/blog/2016/12/practical-guide-to-implement-machine-learning-with-caret-package-in-r-with-practice-problem/

library(tidyverse)
library(caret)
library(randomForest)

rettDmrFull <- read.delim("Rett_sig_individual_smoothed_DMR_methylation.txt")
rettInfo <- read.csv("Rett_sample_info.csv")

# exclude range of columns from 'width' to 'RawDiff', transpose, add diagnosis column
rettDmr <- rettDmrFull %>% 
    as.tibble() %>% 
    select(-(width:RawDiff)) %>%
    unite(seqId1, seqnames, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") %>%
    # transpose: cols to rows
    gather(sampleID, values, -seqId) %>% # cols to rows
    # transpose: rows to cols
    spread(seqId, values) %>% # rows to cols
    # add diagnosis column
    mutate(diagnosis = c("control", "control", "control", 
                         "rett", "rett", "rett", 
                         "rett", "rett", "rett", 
                         "control", "control", "control")) 

rettDmr <- rettDmr %>% select(-sampleID)
View(dmr)



seed <- 9999
set.seed(seed)
trainIndex <- createDataPartition(rettDmr$diagnosis, 
                                  p = 0.8,
                                  list = FALSE )

training <- rettDmr[trainIndex, ]
testing <- rettDmr[-trainIndex, ]


fitControl <- trainControl(method = "repeatedcv", # 10-fold cv
                           number = 10, 
                           repeats = 10)
set.seed(seed)
rf_default <- train( diagnosis ~ ., 
                     data = training, 
                     method = "rf", 
                     trControl = fitControl)
print(rf_default)

# Checking variable importance for random forest
vi <- varImp(object = rf_default)
vi
# Plotting variable importance for random forest
vi_plot <- plot( vi, main = "Random Forest - Variable Importance")
vi_plot

# predictions
predictions <- predict.train( object = rf_default, testing, type = "raw")
table(predictions)

# ----------------------------------------------------------------------

cl <- makePSOCKcluster(4)
registerDoParallel(cl)

# fit a random forest model using ranger 
rf_fit <- train( diagnosis ~ ., 
                 data = training, 
                 method = "ranger" )
stopCluster(cl)

# parallel processing:
# http://topepo.github.io/caret/parallel-processing.html

tRettDmr[1, c(1, 2)]
tRettDmr[1, ]
tRettDmr[ , 1]

#rettDmr $>$
#  rownames_to_column %>% 
#  gather(var, value, -rowname) %>% 
#  spread(rowname, value)



# perform training using random forest
rf_classifier <- randomForest(training, 
                              ntree = 100, 
                              mtry = 2, 
                              importance = TRUE)
rf_classifier

