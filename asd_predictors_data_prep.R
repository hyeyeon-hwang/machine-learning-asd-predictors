# in "asd_predictors_cbdata_prep.R"
# dmrFull, dmrFullCb, sampleInfo  <- list()
# library(tidyverse)
source("asd_predictors_cb_data_prep.R") # dmrFull$cb
source("asd_predictors_mi_data_prep.R") # dmrFull$mi, sampleInfo$mi

# Read data
dmrFull$rett <- read.delim("../data/Individual/Rett_sig_individual_smoothed_DMR_methylation.txt") #, check.names = FALSE)
dmrFull$dup <- read.delim("../data/Individual/Dup15q_sig_individual_smoothed_DMR_methylation.txt")
dmrFull$asd <- read.delim("../data/Individual/ASD_sig_individual_smoothed_DMR_methylation.txt")

dmrFull$plac <- read.delim("../data/Individual/plac_sig_individual_smoothed_DMR_methylation.txt")
dmrFullCb$plac <- read.delim("../data/Consensus_background/background_region_individual_smoothed_methylation.txt")
sampleInfo$plac <- read.csv("../data/Sample_info/sample_info.csv")

#' prepData
#' @description Filter (exclude columns "width" to "RawDiff") and transpose DMR dataset
#' @param dmrFull DMR dataset 
#' @import tidyverse
#' @export prepData
prepData <- function(dmrFull, sampleInfo, dataType = "individual") {
  if(dataType == "individual") {
    dataPrep1 <- dmrFull %>% 
      as.tibble() %>% 
      # instead of select(-(width:RawDiff)), select(-(width:percentDifference)) for plac
      select(-(4:15)) %>%
      unite(seqId1, seqnames, start, sep = ":") %>%
      unite(seqId, seqId1, end, sep = "-") %>%
      gather(sampleID, values, -seqId) %>% # cols to rows
      spread(seqId, values) # rows to cols
    
  }
  
  if(dataType == "cb" | dataType == "mi") {
    dataPrep1 <- dmrFull
    seqId <- row.names(dataPrep1)
    dataPrep1 <- dataPrep1 %>%
      as.tibble() %>%
      add_column(seqId = seqId, .before = 1) %>%
      gather(sampleID, values, -seqId) %>%
      spread(seqId, values)
    
  }
  
  # Add diagnosis column for each sample by matching from sample info files and remove sample ID column
  dataPrep2 <- dataPrep1 %>%  
    add_column(diagnosis = sampleInfo$Diagnosis[match(dataPrep1$sampleID, sampleInfo$Name)], .after = 1)
  return(dataPrep2)
}

dmr <- list()
dmr$rett <- prepData(dmrFull$rett, sampleInfo$rett)
dmr$dup <- prepData(dmrFull$dup, sampleInfo$dup)
dmr$asd <- prepData(dmrFull$asd, sampleInfo$asd)
dmr$plac <- prepData(dmrFull$plac, sampleInfo$plac)
dmr$cb <- prepData(dmrFull$cb, sampleInfo$cb, dataType = "cb")
dmr$mi <- prepData(dmrFull$mi, sampleInfo$mi, dataType = "mi")
dmr$miGrouped <- prepData(dmrFull$mi, sampleInfo$miGrouped, dataType = "mi")
