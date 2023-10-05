############### The following preloading was used for the whole script ###############
# BiocManager::install('curatedMetagenomicData')
# curatedMetagenomicData::sampleMetadata
require(dplyr)
require(ggplot2)
require(curatedMetagenomicData)


############### Download species abundance data from curatedMetagenomicData ############### 
# FengQ <- curatedMetagenomicData("2021-03-31.FengQ_2015.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# GuptaA <- curatedMetagenomicData("2021-03-31.GuptaA_2019.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# HanniganGD <- curatedMetagenomicData("2021-03-31.HanniganGD_2017.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# ThomasAM_2018a <- curatedMetagenomicData("2021-03-31.ThomasAM_2018a.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# ThomasAM_2018b <- curatedMetagenomicData("2021-03-31.ThomasAM_2018b.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# ThomasAM_2019 <- curatedMetagenomicData("2021-03-31.ThomasAM_2019_c.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# VogtmannE <- curatedMetagenomicData("2021-03-31.VogtmannE_2016.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# WirbelJ <- curatedMetagenomicData("2021-03-31.WirbelJ_2018.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# YachidaS <- curatedMetagenomicData("2021-10-14.YachidaS_2019.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# YuJ <- curatedMetagenomicData("2021-03-31.YuJ_2015.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# ZellerG <- curatedMetagenomicData("2021-03-31.ZellerG_2014.relative_abundance", dryrun = FALSE, rownames = 'short')[[1]]
# data <- list(FengQ = FengQ, GuptaA = GuptaA, HanniganGD = HanniganGD, 
#              ThomasAM_2018a = ThomasAM_2018a, ThomasAM_2018b = ThomasAM_2018b, 
#              ThomasAM_2019 = ThomasAM_2019,
#              VogtmannE = VogtmannE, WirbelJ = WirbelJ, YachidaS = YachidaS, YuJ = YuJ, ZellerG = ZellerG)
# 
# saveRDS(data, './case2/data/data_sp.RDS')


############### Data preprocessing ############### 
data <- readRDS('./case2/data/data_sp.RDS')

# Remove low quality cohorts, refer to the <Check data quality> code section for details
data <- data[-c(2,3,6)]

# Get the names of the species detected by all cohorts
common_species <- Reduce(intersect, lapply(data, rownames))

# Combine the relative abundance data of each cohort
relative_abundance <- lapply(data, function(x) x[common_species,]@assays@data$relative_abundance)
relative_abundance <- t(do.call(cbind, relative_abundance))
colnames(relative_abundance) <- common_species

# Log transformation
log_abundance <- log(relative_abundance + 1)

# Rank transformation
N <- nrow(relative_abundance)
P <- ncol(relative_abundance)
ranked_abundance <- (apply(relative_abundance, 1, rank, ties.method = 'min') -1 )/(P-1)
ranked_abundance <- t(ranked_abundance)

# Make sample annotation
cohort <- do.call(c,lapply(data, function(x) x@colData$study_name))
sample_type <- do.call(c,lapply(data, function(x) x@colData$study_condition))
metadata <- data.frame(cohort, sample_type)
rownames(metadata) <- rownames(relative_abundance)

# Keep only control samples and CRC samples 
ranked_abundance <- ranked_abundance %>% subset(. ,sample_type %in% c('control','CRC'))
log_abundance <- log_abundance %>% subset(. ,sample_type %in% c('control','CRC'))
metadata <- metadata[rownames(ranked_abundance),]

# Save the data
saveRDS(list(ranked_abundance, log_abundance, metadata, data),
        './case2/data/processed_data_sp.RDS')


############### Check data quality (not necessary) ############### 
data <- readRDS('./case2/data/data_sp.RDS')
# Get study names
summary_stat <- data.frame(study_names = 
                             names(table(do.call(c,lapply(data, function(x) x@colData$study_name))))
                             )
summary_stat$n_sample <- sapply(data, function(x) ncol(x)) # Get the number of samples
summary_stat$n_species <- sapply(data, function(x) nrow(x)) # Get the number of species
common_species <- Reduce(intersect, lapply(data, rownames)) # Get the names of the species detected by all cohorts
summary_stat$n_comm_species <- length(common_species)
# Get the study condition
summary_stat <- cbind(summary_stat,
                      do.call(bind_rows, lapply(data, function(x) table(x@colData$study_condition))) ) 
# Get the country of each cohort
summary_stat <- cbind(summary_stat,
                      do.call(bind_rows, lapply(data, function(x) table(x@colData$country)))  )
summary_stat[is.na(summary_stat)] <- 0
summary_stat

# Combine relative abundance data of each cohort
relative_abundance <- lapply(data, function(x) x[common_species,]@assays@data$relative_abundance)
relative_abundance <- t(do.call(cbind, relative_abundance))
colnames(relative_abundance) <- common_species

# Check sequencing depth and sparsity
number_reads <- do.call(c, lapply(data, function(x) x@colData$number_reads))
number_sparse <- rowSums(relative_abundance==0) 
cohort <- do.call(c,lapply(data, function(x) x@colData$study_name))

# Visualize sequencing depth and sparsity
ggplot() + geom_boxplot(aes(cohort, number_reads)) + 
  theme(text = element_text(size = 15)) + xlab('Study name')+ ylab('Number of reads')
ggplot() + geom_boxplot(aes(cohort, number_sparse)) +
  theme(text = element_text(size = 15)) + xlab('Number of reads')

# Therefore we exclude the cohorts (GuptaA_2019, HanniganGD_2017, ThomasAM_2019) 
# with very low sequencing depth and high sparsity from our analysis