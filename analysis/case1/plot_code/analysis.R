############### The following preloading was used for Figure 2A (MUST RUN) ############### 
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate/')
load('./case1/r_object/esr1_gene/ESR1.RData')
source('./stablemate.R')

############### Figure 2A ############### 
# As prefiltering of predictors was done before stablemate, we adjust the
# importance score of the pseudo predictor, which was not prefiltered,
# for details on why we did such prefiltering, refer to Supplementary Methods 7.1.6
change_significance <- function(SRST2E_OBJ, n = 100){
  
  # Multiply the importance score of the pseudo predictor by the Lasso selection
  # frequency on the nth most frequently selected predictor of Lasso.
  SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$bt_sel_prob[,1] <- 
    SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$bt_sel_prob[,1] * 
    -sort(-colMeans(SRST2E_OBJ$prediction_ensemble$var_pool))[n]
  
  # Calculate significance for predictivity
  SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$significance <- 
    colMeans(SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$bt_sel_prob > 
               SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$bt_sel_prob[,1])
  
  SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$lt_significance <- 
    colMeans(SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$bt_sel_prob < 
               SRST2E_OBJ$prediction_ensemble$selection_prob$marginal$bt_sel_prob[,1])
  
  SRST2E_OBJ
  
}

# Apply changing of significance
mod <- change_significance(mod)
# Make plot
plot.SRST2E(mod) +  theme(text = element_text(size = 15)) + xlab('Predictivity score') +
  ylab('Stability score')


############### FOR THE FOLLOWING PLOT, CLEAN UP MEMORY, RESTART R AND RERUN PRELOADING ############### 


############### The following preloading was used for Figure 2B,2C,2D,2E S3B, S3C (MUST RUN) ############### 
require(reshape2)
require(dplyr)
require(ggpubr)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate/')
load('./case1/r_object/esr1_pc/ESR1.RData')
source('./stablemate.R')

############### Figure 2B ############### 
plot.SRST2E(mod, 
            label_which = c('PC1','PC2', 'PC3', 'PC4', 'PC5', 'PC6')) +
  theme(text = element_text(size = 15)) + xlab('Predictivity score') +
  ylab('Stability score')


############### Figure 2C ############### 
# Disease type annotation
pam50 <- meta_data$paper_BRCA_Subtype_PAM50
pam50[is.na(pam50)] <- 'Normal'
sample_type <- meta_data$sample_type
sample_type[sample_type == 'Primary Tumor'] <- 'ER+'

# Prepare the data for plotting 
pc_data <- data.frame(type = meta_data$sample_type,
                      pam50 = pam50, pca$x[,1:3],
                      Y = Y, sample_type = sample_type)

# Calculate variance explained by each PC
var_explained <- round(((pca$sdev^2)/sum(pca$sdev^2)) * 100,2)

# Make plot
p1 <- ggplot(data = pc_data) + 
  geom_point(aes(PC1, Y, col = sample_type))+
  geom_smooth(aes(PC1, Y, col = sample_type), method = 'lm') +
  xlab('PC scores') + ylab('ESR1') + 
  scale_color_manual('Sample type', values = c('orange','black')) + 
  theme(text = element_text(size = 15))+
  xlab(paste('PC1', '(', var_explained[1],'% variance explained)'))
p2 <- ggplot(data = pc_data) + 
  geom_point(aes(PC2, Y, col = sample_type))+
  geom_smooth(aes(PC2, Y, col = sample_type), method = 'lm') +
  xlab('PC scores') + ylab('ESR1') + 
  scale_color_manual('Sample type', values = c('orange','black')) + 
  theme(text = element_text(size = 15))+
  xlab(paste('PC2', '(', var_explained[2],'% variance explained)')) + ylab('')
p3 <- ggplot(data = pc_data) + 
  geom_point(aes(PC3, Y, col = sample_type))+
  geom_smooth(aes(PC3, Y, col = sample_type), method = 'lm') +
  xlab('PC scores') + ylab('ESR1') + 
  scale_color_manual('Sample type', values = c('orange','black')) + 
  theme(text = element_text(size = 15))+
  xlab(paste('PC3', '(', var_explained[3],'% variance explained)'))+ ylab('')
ggarrange(p1,p2,p3,common.legend = T, ncol = 3, legend = 'right')


############### Figure S3B ############### 
# Get PC1 loadings and color the first 200 genes
PC1loading <- pca$rotation[,'PC1']
PC1loading <- PC1loading[order(abs(PC1loading), decreasing = T)]
idx <- 1:length(PC1loading)
PC1label <- names(PC1loading)
PC1color <- rep('lightgrey', length(PC1loading))
PC1color[idx <= 200 & PC1loading < 0] <- 'red'
PC1color[idx <= 200 & PC1loading > 0] <- 'blue'

# Make the loading plot for PC1
pload1 <- ggplot() + geom_col(aes(x =  idx,y = PC1loading), col = PC1color) + scale_x_log10() +
  geom_label(aes(x = 2, y =0.03,
                 label = paste(PC1label[1:10], collapse = '\n')),size = 3, color = "red") + 
  geom_label(aes(x = 10, y =0.03,
                 label = paste(PC1label[11:20], collapse = '\n')),size = 3, color = "red") +
  ylim(max(abs(PC1loading)) * c(-1,1)) +
  theme(text = element_text(size =19)) + xlab('Genes ordered by abs PC1 loading') + 
  ylab('PC1 loading') 

# Get PC3 loadings and color the first 200 genes
PC3loading <- pca$rotation[,'PC3']
PC3loading <- PC3loading[order(abs(PC3loading), decreasing = T)]
idx <- 1:length(PC3loading)
PC3label <- names(PC3loading)
PC3color <- rep('lightgrey', length(PC3loading))
PC3color[idx <= 200 & PC3loading < 0] <- 'red'
PC3color[idx <= 200 & PC3loading > 0] <- 'blue'

# Make the loading plot for PC3
pload3 <- ggplot() + geom_col(aes(x =  idx,y = PC3loading), col = PC3color) + scale_x_log10() +
  geom_label(aes(x = 2, y =-0.03,
                 label = paste(PC3label[1:10], collapse = '\n')),size = 3, color = "blue") + 
  geom_label(aes(x = 12, y =-0.03,
                 label = paste(PC3label[11:20], collapse = '\n')),size = 3, color = "blue") +
  ylim(max(abs(PC3loading)) * c(-1,1)) +
  theme(text = element_text(size =19)) + xlab('Genes ordered by abs PC3 loading')+ 
  ylab('PC3 loading')

# Combine the two plots
ggarrange(pload1,pload3,common.legend = T, ncol = 1, legend = 'right')


############### Figure 2D (Running S3B is required) ############### 
# Save top 200 genes for each PC
writeLines(PC1label[1:200], './case1/saved_file/PC1loading.txt')
writeLines(PC3label[1:200], './case1/saved_file/PC3loading.txt')

# The above saved file were used to ran online GO enrichment analysis 
# provided by ShinyGO 0.77 (url:http://bioinformatics.sdstate.edu/go/).
# We downloaded the result returned by ShinyGO, saved as PC1enrichment.csv and
# PC3enrichment.csv for PC1 and PC3 respectively

# Load the GO enrichment result for PC1 and prepare the data for bubble plot
enrichment1 <- read.csv('./case1/downloaded_file/PC1enrichment.csv')
enrichment1 <- enrichment1[1:15,]
enrichment1$Pathway <- factor(enrichment1$Pathway, 
                              levels = enrichment1$Pathway[order(enrichment1$Enrichment.FDR)])

# Load the GO enrichment result for PC1 and prepare the data for bubble plot
enrichment3 <- read.csv('./case1/downloaded_file/PC3enrichment.csv')
enrichment3 <- enrichment3[1:15,]
enrichment3$Pathway <- factor(enrichment3$Pathway,
                              levels = enrichment3$Pathway[order(enrichment3$Enrichment.FDR)])

# Make the loading plot for PC1
penrich1 <- ggplot(data = enrichment1) + geom_point(aes(x = -log10(Enrichment.FDR),
                                                        y = Pathway, size = Fold.Enrichment,
                                                        col = -log10(Enrichment.FDR))) +
  scale_color_gradient('Negative log10 FDR', low = 'blue',high = 'red') +
  scale_size_continuous('Fold enrichment')+
  theme(text = element_text(size = 16)) + 
  xlab('Negative log10 FDR') + ylab('Pathway') +
  guides(color = 'none') 

# Make the loading plot for PC3
penrich3 <- ggplot(data = enrichment3) + geom_point(aes(x = -log10(Enrichment.FDR),
                                                        y = Pathway, size = Fold.Enrichment,
                                                        col = -log10(Enrichment.FDR))) +
  scale_color_gradient('Negative log10 FDR', low = 'blue',high = 'red') +
  scale_size_continuous('Fold enrichment')+
  theme(text = element_text(size = 16)) + 
  xlab('Negative log10 FDR') + ylab('Pathway') +
  guides(color = 'none') 

# Combine the two plots
ggarrange(penrich1,penrich3,common.legend = T, ncol = 1, legend = 'right',align = 'v')


############### Figure 2E (Running S3B is required) ############### 
# Function for aggregating gene expression, refer to Section 4.3 for details
npca <- function(x){
  pca <- prcomp(x, center = F, scale = F)
  l <- pca$rotation[,1]
  l <- sign(l[which.max(abs(l))]) * l
  l <- l - min(l)
  l <- l/sqrt(sum(l^2))
  as.numeric(x %*% l)
}

# Load the METABRIC data
ma_data <- read.delim('./case1/data/cBioProtal/brca_metabric/data_mrna_agilent_microarray_zscores_ref_all_samples.txt')
row.names(ma_data) <- ma_data$Hugo_Symbol
ma_data <- ma_data[,-2:-1]

# Load the METABRIC annotation
ma_metadata <- read.delim('./case1/data/cBioProtal/brca_metabric/data_clinical_patient.txt')
rownames(ma_metadata) <- ma_metadata$X.Patient.Identifier

# Filter the METABRIC data to keep only the ER+ samples
ma_pam50 <- ma_metadata[gsub('\\.','-',colnames(ma_data)),]$Pam50...Claudin.low.subtype
ma_data <- ma_data[,ma_pam50 %in% c('LumA','LumB')]
ma_data[is.na(ma_data)] <- 0
ma_pam50 <- ma_pam50[ma_pam50 %in% c('LumA','LumB')]

# Filter the METABRIC data to keep only the genes shared with the top 200 genes of PC3 
ma_common_gene <- intersect(rownames(ma_data), PC3label[1:200])
ma_scores <- t(ma_data[ma_common_gene,]) %>% npca # Aggregate gene expression
ma_esr1 <- as.numeric(ma_data['ESR1',]) # Extract ESR1 

# Load the Gtex data
gtex_data <- read.delim('./case1/data/GTEx/gene_tpm_2017-06-05_v8_breast_mammary_tissue.gct',
                        skip = 2)
gtex_gname <- gtex_data$Description
gtex_data <- gtex_data[,-3:-1]
gtex_data <- log(gtex_data + 1) # log transform gtex data

# Clean up the names of the Gtex samples 
gtex_sname <- colnames(gtex_data)
gtex_sname <- substr(gtex_sname, 1,  nchar(gtex_sname)-14)
gtex_sname <- gsub('\\.','-',gtex_sname)

# Load the Gtex annotation
gtex_meta_data <- read.delim('./case1/data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
rownames(gtex_meta_data) <- gtex_meta_data$SUBJID

# Filter the Gtex data to keep only the genes shared with the top 200 genes of PC3 
gtex_common_gene <- intersect(gtex_gname, PC3label[1:200])
gtex_scores <- t(gtex_data[gtex_gname %in% gtex_common_gene,]) %>% npca # Aggregate gene expression
gtex_esr1 <-   as.numeric(gtex_data[gtex_gname == 'ESR1',]) # Extract ESR1 

# Prepare the data for plotting 
ex_data <- data.frame(esr1 = c(ma_esr1, gtex_esr1),
                      scores = c(ma_scores, gtex_scores),
                      type = c(rep('Cancer', length(ma_scores)), rep('Normal', length(gtex_scores))),
                      sample_type = c(rep('ER+', length(ma_scores)), rep('Normal', length(gtex_scores))),
                      pam50 = c(ma_pam50, rep('Normal', length(gtex_scores))),
                      data = c(rep('METABRIC', length(ma_scores)), rep('GTEx', length(gtex_scores))))

# Make the plot
ggplot(data = ex_data) + geom_point(aes(scores, esr1, col = sample_type), show.legend = T) +
  geom_smooth(aes(scores, esr1, col = sample_type), method = 'lm')+
  facet_wrap(~data, scales = 'free', ncol = 1)+
  xlab('Aggregated expression of the top 200 PC3 genes') + ylab('ESR1') + 
  scale_color_manual('Sample type', values = c('orange','black')) + 
  theme(text = element_text(size = 15), axis.title.x = element_text(hjust = 0.2))  


############### Figure S3C ############### 
# Function for aggregating gene expression, refer to Section 4.3 for details
npca <- function(x){
  pca <- prcomp(x, center = F, scale = F)
  l <- pca$rotation[,1]
  l <- sign(l[which.max(abs(l))]) * l
  l <- l - min(l)
  l <- l/sqrt(sum(l^2))
  as.numeric(x %*% l)
}

# Load the basal gene list provided by Li et al. (2021)
basal_gene <- readLines('./case1/downloaded_file/basal_gene.txt')
basal_gene <- intersect(basal_gene, colnames(X)) # Only consider the genes that are also measured by TCGA

# Disease type annotation
pam50 <- meta_data$paper_BRCA_Subtype_PAM50
pam50[is.na(pam50)] <- 'Normal'
sample_type <- meta_data$sample_type
sample_type[sample_type == 'Primary Tumor'] <- 'ER+'
basal_score <- X[,basal_gene] %>% npca

# Plot the aggregated expression of basal genes against ESR1
p1 <- ggplot() + geom_point(aes(basal_score, Y, col = sample_type)) + 
  geom_smooth(aes(basal_score, Y, col = sample_type), method = 'lm')+
  scale_color_manual(values = c('orange','black')) + xlab('') + ylab('ESR1') +
  theme(text = element_text(size = 18)) + guides(col = F)

# Plot the aggregated expression of basal genes against PC3
p2 <- ggplot() + geom_point(aes(basal_score, pca$x[,3], col = sample_type)) + 
  geom_smooth(aes(basal_score, pca$x[,3], col = sample_type), method = 'lm')+
  scale_color_manual('Sample type', values = c('orange','black'))  + xlab('') + ylab('PC3') +
  theme(text = element_text(size = 18), axis.title.x = element_text(hjust = -0.3))

# Combine the two plots
p <- ggarrange(p1,p2, common.legend = T, nrow = 1, legend = 'right')
annotate_figure(p, bottom = text_grob('Aggregated expression of basal BC markers', 
                                      size = 18,vjust = -1, hjust = 0.6))
