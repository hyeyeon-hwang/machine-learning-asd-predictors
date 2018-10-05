library(tidyverse)
library(caret)
library(randomForest)


rettDmrFull <- read.delim("Rett_sig_individual_smoothed_DMR_methylation.txt")
rettInfo <- read.csv("Rett_sample_info.csv")



# exclude range of columns from 'width' to 'RawDiff'
rettDmr <- rettDmrFull %>% as.tibble() %>% select(-(width:RawDiff)) 
rettDmr

rettDmr <- unite(rettDmr, 
                 seqId1, 
                 seqnames, start, 
                 sep = ":")

rettDmr <- unite(rettDmr, 
                 seqId, 
                 seqId1, end, 
                 sep = "-")

# as.tibble( t(rettDmr) )
tRettDmr <- rettDmr %>% t() %>% as.tibble()

tRettDmr$diagnosis <- c("diagnosis", "control", "control", "control", 
                        "rett", "rett", "rett", 
                        "rett", "rett", "rett", 
                        "control", "control", "control")

set.seed(5)
trainIndex <- createDataPartition(tRettDmr$diagnosis, 
                                  p = 0.8,
                                  list = FALSE )

training <- tRettDmr[trainIndex, ]
testing <- tRettDmr[-trainIndex, ]

cl <- makePSOCKcluster(4)
registerDoParallel(cl)
# fit a random forest model using ranger 
rf_fit <- train( as.factor(diagnosis) ~ ., 
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

