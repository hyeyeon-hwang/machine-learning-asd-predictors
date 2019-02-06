library(tidyverse)

dmrFullCb <- list()
dmrFullCb$rett <- read.delim("../data/Consensus_background/Rett_consensus_background_individual_smoothed_DMR_methylation.txt", check.names = FALSE)
dmrFullCb$dup <- read.delim("../data/Consensus_background/Dup_consensus_background_individual_smoothed_DMR_methylation.txt")
dmrFullCb$asd <- read.delim("../data/Consensus_background/ASD_consensus_background_individual_smoothed_DMR_methylation.txt")

sampleInfo <- list()
sampleInfo$rett <- read.csv("../data/Sample_info/Rett_sample_info.csv") %>% add_column(batch = 1)
sampleInfo$dup <- read.csv("../data/Sample_info/Dup15q_sample_info.csv") %>% add_column(batch = 2)
sampleInfo$asd <- read.csv("../data/Sample_info/ASD_sample_info.csv") %>% add_column(batch = 3)

sampleInfo$cb <- tibble(Name = c(as.character(sampleInfo$rett$Name), 
                                         as.character(sampleInfo$dup$Name),
                                         as.character(sampleInfo$asd$Name)), 
                            Diagnosis = c(as.character(sampleInfo$rett$Diagnosis), 
                                          as.character(sampleInfo$dup$Diagnosis), 
                                          as.character(sampleInfo$asd$Diagnosis)), 
                            batch = c(as.character(sampleInfo$rett$batch), 
                                      as.character(sampleInfo$dup$batch), 
                                      as.character(sampleInfo$asd$batch)))

preCombatPrep <- function(dmrFullCb) {
  data <- dmrFullCb %>%
    as.tibble() %>% 
    unite(seqId1, chr, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") 
  return(data)
}

dmrPreCombat <- list()
dmrPreCombat$rett <- preCombatPrep(dmrFullCb$rett)
dmrPreCombat$dup <- preCombatPrep(dmrFullCb$dup)
# remove repeated samples: JLKD063 = 1136 , JLKD066 = 1406, JLKD067 = 1711
dmrPreCombat$asd <- preCombatPrep(dmrFullCb$asd) %>% select(-c("JLKD063", "JLKD066", "JLKD067"))

# joinedCb for combat  
dmrPreCombat$joinedCb <- dmrPreCombat$rett %>%
  full_join(dmrPreCombat$dup, by = "seqId") %>%
  full_join(dmrPreCombat$asd, by = "seqId") %>%
  drop_na() %>%
  gather(sampleID, values, -seqId) %>% # transpose: cols to rows
  spread(seqId, values)# transpose: rows to cols
dmrPreCombat$joinedCb <- dmrPreCombat$joinedCb %>%
  add_column(diagnosis = as.factor(sampleInfo$cb$Diagnosis[match(dmrPreCombat$joinedCb$sampleID, sampleInfo$cb$Name)]), .after = 1) %>%
  add_column(batch = as.numeric(sampleInfo$cb$batch[match(dmrPreCombat$joinedCb$sampleID, sampleInfo$cb$Name)]), .after = 2) %>%
  add_column(sample = as.integer(1:nrow(dmrPreCombat$joinedCb)), .before = 1)
groupedDiagnosis <- as.character(dmrPreCombat$joinedCb$diagnosis)
groupedDiagnosis[which(groupedDiagnosis != "Control")] <- "Positive"

dmrPreCombat$joinedCb <- dmrPreCombat$joinedCb %>% 
  add_column(groupedDiagnosis = as.factor(groupedDiagnosis), .after = 1)

modCombatData <- dmrPreCombat$joinedCb[, 1:5]
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
postCombatJoinedCb = ComBat(dat = datCombat, batch = batch, mod = modCombat, par.prior = TRUE, prior.plots = FALSE)
#View(postCombatJoinedCb)
#write.table(postCombatJoinedCb, "postCombatJoinedCb.txt", sep = "\t")

dmrFull <- list()
dmrFull$cb <- postCombatJoinedCb
