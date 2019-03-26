library(pheatmap)

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")
example_file <- system.file ("extra/TagSeqExample.tab", package="DESeq")
data <- read.delim(example_file, header=T, row.names="gene")
data_subset <- as.matrix(data[rowSums(data)>50000,])
data_subset[1:4, 1:2]

pheatmap(data_subset[1:4, 1:2])

rfModel$mi3 # 115  0.6250000  0.4137285
rfModel$mi3$finalModel
rfModel$mi3$pred
cmMi3Perc <- confusionMatrix.train(rfModel$mi3)
cmMi3Cnt <- confusionMatrix.train(rfModel$mi3, norm = "none")

library(gplots)
mergedBrainAcc <- data.frame(
  "Class Accuracy" = c(10/(10 + 6),
                       14/(5+14),
                       2/(2+3),
                       3/(3+3)),
  "Class" = c("ASD","Control","Dup15q","Rett"))
heatmap.2(cbind(mergedBrainAcc$Class.Accuracy, mergedBrainAcc$Class.Accuracy),
          trace = "n",
          Colv = NA,
          dendrogram = "none",
          labCol = "",
          labRow = mergedBrainAcc$Class,
          cexRow = 1,
          main = "Merged Brain Class Accuracies",
          xlab = "Class Accuracy",
          ylab = "Class")
  # c(
  # "ASD" = 10/(10 + 6),
  # "Control" = 14/(5+14),
  # "Dup15q" = 2/(2+3),
  # "Rett" = 3/(3+3)
  # )

mBrainAcc <- matrix(c(10/(10 + 6),
                       14/(5+14),
                       2/(2+3),
                       3/(3+3))) %>%
  t()

colors <- c("floralwhite","lightsalmon","indianred1","firebrick3","firebrick","firebrick4")
colnames(mBrainAcc) <- c("ASD","Control","Dup15q","Rett")

pdf("4classModel.pdf", height = 4, width = 11, onefile = FALSE)
pheatmap(mBrainAcc, cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(colors)(50),
         main = "4-class Random Forest Model")
dev.off()

