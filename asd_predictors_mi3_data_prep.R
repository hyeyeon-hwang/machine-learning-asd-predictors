library(tidyverse)

dmrFullMi3 <- list()
dmrFullMi3$rett <- read.delim("../data/Merged_individual/fromMerged3_asd_dup_rett/rett_merged_individual_dmr.txt") #, check.names = FALSE)
dmrFullMi3$dup <- read.delim("../data/Merged_individual/fromMerged3_asd_dup_rett/dup_merged_individual_dmr.txt")
dmrFullMi3$asd <- read.delim("../data/Merged_individual/fromMerged3_asd_dup_rett/asd_merged_individual_dmr.txt")

# sampleInfo <- list(), already in ...cb_data_prep.R
sampleInfo$rett <- read.csv("../data/Sample_info/Rett_sample_info.csv") %>% add_column(batch = 1)
sampleInfo$dup <- read.csv("../data/Sample_info/Dup15q_sample_info.csv") %>% add_column(batch = 2)
sampleInfo$asd <- read.csv("../data/Sample_info/ASD_sample_info.csv") %>% add_column(batch = 3)

# put X's in front of sample name in rett 
sampleInfo$rett$Name <- paste("X", sampleInfo$rett$Name, sep = "")

sampleInfo$mi3 <- tibble(Name = c(as.character(sampleInfo$rett$Name), 
                                 as.character(sampleInfo$dup$Name),
                                 as.character(sampleInfo$asd$Name)), 
                        Diagnosis = c(as.character(sampleInfo$rett$Diagnosis), 
                                      as.character(sampleInfo$dup$Diagnosis), 
                                      as.character(sampleInfo$asd$Diagnosis)), 
                        batch = c(as.character(sampleInfo$rett$batch), 
                                  as.character(sampleInfo$dup$batch), 
                                  as.character(sampleInfo$asd$batch)))

preCombatPrep <- function(dmrFullMi3) {
  data <- dmrFullMi3 %>%
    as.tibble() %>% 
    unite(seqId1, chr, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") 
  return(data)
}

dmrPreCombat <- list()
dmrPreCombat$rett <- preCombatPrep(dmrFullMi3$rett)
dmrPreCombat$dup <- preCombatPrep(dmrFullMi3$dup)
# remove repeated samples: JLKD063 = 1136 , JLKD066 = 1406, JLKD067 = 1711
dmrPreCombat$asd <- preCombatPrep(dmrFullMi3$asd) %>% select(-c("JLKD063", "JLKD066", "JLKD067"))

# joinedCb for combat  
dmrPreCombat$joinedMi3 <- dmrPreCombat$rett %>%
  full_join(dmrPreCombat$dup, by = "seqId") %>%
  full_join(dmrPreCombat$asd, by = "seqId") %>%
  drop_na() %>%
  gather(sampleID, values, -seqId) %>% # transpose: cols to rows
  spread(seqId, values)# transpose: rows to cols
dmrPreCombat$joinedMi3 <- dmrPreCombat$joinedMi3 %>%
  add_column(diagnosis = as.factor(sampleInfo$mi3$Diagnosis[match(dmrPreCombat$joinedMi3$sampleID, sampleInfo$mi3$Name)]), .after = 1) %>%
  add_column(batch = as.numeric(sampleInfo$mi3$batch[match(dmrPreCombat$joinedMi3$sampleID, sampleInfo$mi3$Name)]), .after = 2) %>%
  add_column(sample = as.integer(1:nrow(dmrPreCombat$joinedMi3)), .before = 1)
groupedDiagnosis <- as.character(dmrPreCombat$joinedMi3$diagnosis)
groupedDiagnosis[which(groupedDiagnosis != "Control")] <- "Positive"

dmrPreCombat$joinedMi3 <- dmrPreCombat$joinedMi3 %>% 
  add_column(groupedDiagnosis = as.factor(groupedDiagnosis), .after = 1)

modCombatData <- dmrPreCombat$joinedMi3[, 1:5]
order <- c("sampleID", "sample", "diagnosis", "batch", "groupedDiagnosis")
modCombatData <- modCombatData[, order] %>% 
  as.data.frame()
row.names(modCombatData) <- modCombatData$sampleID
modCombatData <- modCombatData[,2:5]

# info_joinedCB needs to be in current format for modcombat
batch <- modCombatData$batch 

datCombat <- dmrPreCombat$rett %>%
  full_join(dmrPreCombat$dup, by = "seqId") %>%
  full_join(dmrPreCombat$asd, by = "seqId") %>%
  drop_na() %>% as.data.frame()
row.names(datCombat) <- datCombat$seqId
datCombat <- datCombat[,-1] %>% as.matrix.data.frame()

# combat modcombat
modCombat = model.matrix(~1, data = modCombatData) #joined_info)

library(sva)
# combat function
postCombatJoinedMi3 = ComBat(dat = datCombat, batch = batch, mod = modCombat, par.prior = TRUE, prior.plots = FALSE)
#View(postCombatJoinedCb)
#write.table(postCombatJoinedCb, "postCombatJoinedCb.txt", sep = "\t")

# dmrFull <- list(), already in ...cb_data_prep.R
dmrFull$mi3 <- postCombatJoinedMi3
sampleInfo$mi3Grouped <- sampleInfo$mi3
sampleInfo$mi3Grouped$Diagnosis[ which(sampleInfo$mi3Grouped$Diagnosis != "Control")] <- "Positive"

  
  