############### The following preloading was used for the whole script ###############
# BiocManager::install('curatedMetagenomicData')
# curatedMetagenomicData::sampleMetadata
require(dplyr)
require(curatedMetagenomicData)
  
############### Download pathway abundance data from curatedMetagenomicData ############### 
# FengQ <- curatedMetagenomicData("2021-03-31.FengQ_2015.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(FengQ)) # For shorter pathway names
# FengQ <- FengQ[-idx,]
# GuptaA <- curatedMetagenomicData("2021-03-31.GuptaA_2019.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(GuptaA))
# GuptaA <- GuptaA[-idx,]
# HanniganGD <- curatedMetagenomicData("2021-03-31.HanniganGD_2017.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(HanniganGD))
# HanniganGD <- HanniganGD[-idx,]
# ThomasAM_2018a <- curatedMetagenomicData("2021-03-31.ThomasAM_2018a.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(ThomasAM_2018a))
# ThomasAM_2018a <- ThomasAM_2018a[-idx,]
# ThomasAM_2018b <- curatedMetagenomicData("2021-03-31.ThomasAM_2018b.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(ThomasAM_2018b))
# ThomasAM_2018b <- ThomasAM_2018b[-idx,]
# ThomasAM_2019 <- curatedMetagenomicData("2021-03-31.ThomasAM_2019_c.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(ThomasAM_2019))
# ThomasAM_2019 <- ThomasAM_2019[-idx,]
# VogtmannE <- curatedMetagenomicData("2021-03-31.VogtmannE_2016.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(VogtmannE))
# VogtmannE <- VogtmannE[-idx,]
# WirbelJ <- curatedMetagenomicData("2021-03-31.WirbelJ_2018.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(WirbelJ))
# WirbelJ <- WirbelJ[-idx,]
# YachidaS <- curatedMetagenomicData("2021-10-14.YachidaS_2019.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(YachidaS))
# YachidaS <- YachidaS[-idx,]
# YuJ <- curatedMetagenomicData("2021-03-31.YuJ_2015.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(YuJ))
# YuJ <- YuJ[-idx,]
# ZellerG <- curatedMetagenomicData("2021-03-31.ZellerG_2014.pathway_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# idx <- grep('\\|',rownames(ZellerG))
# ZellerG <- ZellerG[-idx,]
# data <- list(FengQ = FengQ, GuptaA = GuptaA, HanniganGD = HanniganGD, 
#              ThomasAM_2018a = ThomasAM_2018a, ThomasAM_2018b = ThomasAM_2018b, 
#              ThomasAM_2019 = ThomasAM_2019,
#              VogtmannE = VogtmannE, WirbelJ = WirbelJ, YachidaS = YachidaS, YuJ = YuJ, ZellerG = ZellerG)

# saveRDS(data, './case2/data/data_pw.RDS')

            
############### Data preprocessing ############### 
data <- readRDS('./case2/data/data_pw.RDS')

# Remove low quality cohorts, refer to the <Check data quality> code section in <preprocess_sp.R> for details
data <- data[-c(2,3,6)]

# Get the names of the species detected by all cohorts
common_species <- Reduce(intersect, lapply(data, rownames))

# Combine the relative abundance data of each cohort
pathway_abundance <- lapply(data, function(x) x[common_species,]@assays@data$pathway_abundance)
pathway_abundance <- t(do.call(cbind, pathway_abundance))
colnames(pathway_abundance) <- common_species

# Log transformation
log_abundance <- log(pathway_abundance + 1)

# Rank transformation
N <- nrow(pathway_abundance)
P <- ncol(pathway_abundance)
ranked_abundance <- (apply(pathway_abundance, 1, rank, ties.method = 'min') -1 )/(P-1)
ranked_abundance <- t(ranked_abundance)

# Make sample annotation
cohort <- do.call(c,lapply(data, function(x) x@colData$study_name))
sample_type <- do.call(c,lapply(data, function(x) x@colData$study_condition))
metadata <- data.frame(cohort, sample_type)
rownames(metadata) <- rownames(pathway_abundance)

# Keep only control samples and CRC samples 
ranked_abundance <- ranked_abundance %>% subset(. ,sample_type %in% c('control','CRC'))
log_abundance <- log_abundance %>% subset(. ,sample_type %in% c('control','CRC'))
metadata <- metadata[rownames(ranked_abundance),]

# Save the data
saveRDS(list(ranked_abundance, log_abundance, metadata, data),
        './case2/data/processed_data_pw.RDS')