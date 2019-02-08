library(tidyverse)

dmrFullMi4 <- list()
dmrFullMi4$rett <- read.delim("../data/Merged_individual/fromMerged4_asd_dup_rett_plac/rett_merged4_individual_dmr.txt") #, check.names = FALSE)
dmrFullMi4$dup <- read.delim("../data/Merged_individual/fromMerged4_asd_dup_rett_plac/dup_merged4_individual_dmr.txt")
dmrFullMi4$asd <- read.delim("../data/Merged_individual/fromMerged4_asd_dup_rett_plac/asd_merged4_individual_dmr.txt")
dmrFullMi4$plac <- read.delim("../data/Merged_individual/fromMerged4_asd_dup_rett_plac/plac_merged4_individual_dmr.txt")

# sampleInfo <- list(), already in ...cb_data_prep.R
# sampleInfo$rett <- read.csv("../data/Sample_info/Rett_sample_info.csv") %>% add_column(batch = 1)
# sampleInfo$dup <- read.csv("../data/Sample_info/Dup15q_sample_info.csv") %>% add_column(batch = 2)
# sampleInfo$asd <- read.csv("../data/Sample_info/ASD_sample_info.csv") %>% add_column(batch = 3)
sampleInfo$plac <- read.csv("../data/Sample_info/sample_info.csv") %>% add_column(batch = 4)

# put X's in front of sample name in rett 
#sampleInfo$rett$Name <- paste("X", sampleInfo$rett$Name, sep = "")

sampleInfo$mi4 <- tibble(Name = c(as.character(sampleInfo$rett$Name), 
                                  as.character(sampleInfo$dup$Name),
                                  as.character(sampleInfo$asd$Name),
                                  as.character(sampleInfo$plac$Name)), 
                         Diagnosis = c(as.character(sampleInfo$rett$Diagnosis), 
                                       as.character(sampleInfo$dup$Diagnosis), 
                                       as.character(sampleInfo$asd$Diagnosis),
                                       as.character(sampleInfo$plac$Diagnosis)), 
                         batch = c(as.character(sampleInfo$rett$batch), 
                                   as.character(sampleInfo$dup$batch), 
                                   as.character(sampleInfo$asd$batch),
                                   as.character(sampleInfo$plac$batch)))

preCombatPrep <- function(dmrFullMi4) {
  data <- dmrFullMi4 %>%
    as.tibble() %>% 
    unite(seqId1, chr, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") 
  return(data)
}

dmrPreCombat <- list()
dmrPreCombat$rett <- preCombatPrep(dmrFullMi4$rett)
dmrPreCombat$dup <- preCombatPrep(dmrFullMi4$dup)
# remove repeated samples: JLKD063 = 1136 , JLKD066 = 1406, JLKD067 = 1711
dmrPreCombat$asd <- preCombatPrep(dmrFullMi4$asd) %>% select(-c("JLKD063", "JLKD066", "JLKD067"))
dmrPreCombat$plac <- preCombatPrep(dmrFullMi4$plac)

# joinedCb for combat  
dmrPreCombat$joinedMi4 <- dmrPreCombat$rett %>%
  full_join(dmrPreCombat$dup, by = "seqId") %>%
  full_join(dmrPreCombat$asd, by = "seqId") %>%
  full_join(dmrPreCombat$plac, by = "seqId") %>%
  drop_na() %>%
  gather(sampleID, values, -seqId) %>% # transpose: cols to rows
  spread(seqId, values)# transpose: rows to cols
dmrPreCombat$joinedMi4 <- dmrPreCombat$joinedMi4 %>%
  add_column(diagnosis = as.factor(sampleInfo$mi4$Diagnosis[match(dmrPreCombat$joinedMi4$sampleID, sampleInfo$mi4$Name)]), .after = 1) %>%
  add_column(batch = as.numeric(sampleInfo$mi4$batch[match(dmrPreCombat$joinedMi4$sampleID, sampleInfo$mi4$Name)]), .after = 2) %>%
  add_column(sample = as.integer(1:nrow(dmrPreCombat$joinedMi4)), .before = 1)
groupedDiagnosis <- as.character(dmrPreCombat$joinedMi4$diagnosis)
groupedDiagnosis[which(groupedDiagnosis != "Control")] <- "Positive"

dmrPreCombat$joinedMi4 <- dmrPreCombat$joinedMi4 %>% 
  add_column(groupedDiagnosis = as.factor(groupedDiagnosis), .after = 1)

modCombatData <- dmrPreCombat$joinedMi4[, 1:5]
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
  full_join(dmrPreCombat$plac, by = "seqId") %>%
  drop_na() %>% as.data.frame()
row.names(datCombat) <- datCombat$seqId
datCombat <- datCombat[,-1] %>% as.matrix.data.frame()

# combat modcombat
modCombat = model.matrix(~1, data = modCombatData) #joined_info)

library(sva)
# combat function
postCombatJoinedMi4 = ComBat(dat = datCombat, batch = batch, mod = modCombat, par.prior = TRUE, prior.plots = FALSE)
#View(postCombatJoinedCb)
#write.table(postCombatJoinedCb, "postCombatJoinedCb.txt", sep = "\t")

# dmrFull <- list(), already in ...cb_data_prep.R
dmrFull$mi4 <- postCombatJoinedMi4
sampleInfo$mi4Grouped <- sampleInfo$mi4
sampleInfo$mi4Grouped$Diagnosis[ which(sampleInfo$mi4Grouped$Diagnosis != "Control")] <- "Positive"


