############### The following preloading was used for Figure 4A; S8A,B (MUST RUN) ############### 
require(Sincast)
require(destiny)
require(Seurat)
require(plotly)
require(pheatmap)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate/')
source('./stablemate.R')

# Load the reference data, the query data for Sincast projection
reference <- readRDS('./case3/data/atlas_data/atlas.rds') # Load the atlas data
# Load the Seurat object storing the log-normalized query data and the query annotation
seurat_obj <- readRDS('./case3/data/processed_query_data/seurat_data.rds')
# Load the Sincast imputed query data
impdata <- readRDS('./case3/data/processed_query_data/impdata.rds')
# Load the diffusion pseudo time object
DPT <- readRDS('./case3/data/processed_query_data/dpt.rds')

# The diffusion pseudo time starting from the 366th cell (the cell with the smallest 
# DC1 value) is chosen as the response for stablemate analysis.
DPT366 <- DPT$DPT366

# Get cluster id of cells
location <- seurat_obj$Cluster_2d
# Rename cluster id to cell locations, which will be used to annotate query cells on the reference PCA
location[location ==7] <- 'Periphery'
location[location =='8'] <- 'Core'
# # Get the inferred branching of cells from diffusion pseudo time analysis 
# branch <- DPT@branch[,1]
# # Treat branching of cells as cell locations as well
# location[branch == 3] <- 'Tip2'
# location[branch == 2] <- 'Tip1'


############### Figure 4A (Sincast projection, can be skipped if Sincast is not installed) ############### 
# Save the imputed query data as a sce object
query <- createSce(data = impdata, colData = seurat_obj@meta.data)
query$location <- location

# Weight reference genes by their ability in classifying reference celltypes
reference <- featureWeighting(reference, clusterid = 'celltype', cut = F)

# Filter both the reference and the query genes down to the genes shared by both the data
c(reference, query) %<-% filterData(reference, query)

# Create the reference color scheme
reference_cols <- c("#081d58","#225ea8","#1d91c0","#253494",
                    "#7fcdbb","#c7e9b4","#edf8b1","#41b6c4","lightgrey",
                    "lightgrey","lightgrey","lightgrey","lightgrey","lightgrey",
                    "lightgrey","lightgrey","lightgrey","lightgrey","lightgrey",
                    "lightgrey","lightgrey","lightgrey","lightgrey")
names(reference_cols) <-  c("kupffer cell","microglia","macrophage",
                            "monocyte","CD141+ dendritic cell","CD1c+ dendritic cell",
                            "plasmacytoid dendritic cell","dendritic cell",
                            "common myeloid progenitor","common lymphoid progenitor",
                            "granulocyte monocyte progenitor",
                            "hematopoietic multipotent progenitor","neutrophil",
                            "granulocyte","myelocyte","metamyelocyte","promyelocyte",
                            "erythrocyte","erythroblast","proerythroblast",
                            "endothelial progenitor","hemogenic endothelium",
                            "hemangioblast")

# Create the query color scheme
query_cols <-  c('lavenderblush4', 'black')
names(query_cols)  <- c('Core', 'Periphery')

# Run PCA on the reference data to create an atlas
reference <- makeAtlas(reference = reference, col.by = 'celltype', colors = reference_cols, vis.atlas = T)

# Project the query onto the atlas
query <- project(reference, query, assay = 'data')

# Diffusion reconstruction of the projection, which is what is illustrated in Figure 4A1 where cells are annotated by cell location
DiffusionReconstruct(reference, query,  colReference.by = 'celltype', a = 10,
                     referenceColors = reference_cols, colQuery.by = 'location',
                     queryColors = query_cols)

# Diffusion reconstruction of the projection, which is what is illustrated in Figure 4A2 where cells are annotated by diffusion pseudo time
DiffusionReconstruct(reference, query,  colReference.by = 'celltype', a = 10,
                     referenceColors = reference_cols, colQuery.by = DPT366)


############### Figure S8A (Running the previous code section is required) ############### 
# Capybara cell identity profiling
query <- SincastCapybara(reference, query, clusterid = 'celltype', w = 'HD_mean')

# Heatmap visualization of Capybara cell scores
ano_col <- data.frame(Location = location, DPT = DPT366) # Create column (cell) annotation
rownames(ano_col) <- colnames(query)
# Make the heatmap, order cells by their diffusion pseudo time
pheatmap::pheatmap(t(query@metadata$Capybara[order(DPT366),]),
                   cluster_cols = F, show_colnames = F, annotation_col = ano_col,
                   color = viridis::viridis(20),
                   annotation_colors = list(Losition = c(Core = 'lavenderblush4',
                                                         Periphery = 'black')),
                   treeheight_row = 20, fontsize = 12)


############### Figure S8B ############### 
# Diffusion map visualization, which is what is illustrated in Figure 4B1 where cells are annotated by cell location
ggplot(data = data.frame(attr(DPT@dm,'eigenvectors'))) + 
  geom_point(aes(DC1,DC2,color = location)) + 
  scale_color_manual('Location', values = c('lavenderblush4','black')) +
  theme(legend.position = 'bottom')

# Diffusion map visualization, which is what is illustrated in Figure 4B2 where cells are annotated by diffusion pseudotime
ggplot(data = data.frame(attr(DPT@dm,'eigenvectors'))) + 
  geom_point(aes(DC1,DC2,color = DPT366)) + scale_color_viridis_b('DPT')+
  theme(legend.position = 'bottom')

# The diffusion pseudo time starting from the 366th cell (the cell with the smallest 
# DC1 value) is chosen for analysis.


############### FOR THE FOLLOWING PLOT, CLEAN UP MEMORY, RESTART R AND RERUN PRELOADING ############### 


############### The following preloading was used for Figure 4B-F; S8C-F, S9 (MUST RUN) ############### 
require(destiny)
require(Seurat)
require(pheatmap)
require(ggpubr)
require(mgcv)
require(ggplot2)
require(RColorBrewer)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate/')
source('./stablemate.R')

# Load the Seurat object storing the log-normalized query data and the query annotation
seurat_obj <- readRDS('./case3/data/processed_query_data/seurat_data.rds')
# Load the diffusion pseudo time object
DPT <- readRDS('./case3/data/processed_query_data/dpt.rds')
# The diffusion pseudo time starting from the 366th cell (the cell with the smallest 
# DC1 value) is chosen as the response for stablemate analysis.
DPT366 <- DPT$DPT366

# Get cluster id of cells
location <- seurat_obj$Cluster_2d
# Rename cluster id to cell locations, which will be used to annotate query cells on the reference PCA
location[location ==7] <- 'Periphery'
location[location =='8'] <- 'Core'
# # Get the inferred branching of cells from diffusion pseudo time analysis 
# branch <- DPT@branch[,1]
# # Treat branching of cells as cell locations as well
# location[branch == 3] <- 'Tip2'
# location[branch == 2] <- 'Tip1'

# Load the stablemate result for genes that are predictive in the core
load('./case3/r_object/dpt_gene/core.RData')
mod_core <- mod
# Load the stablemate result for genes that are predictive in the periphery
load('./case3/r_object/dpt_gene/periphery.RData')
mod_periphery <- mod
# Load the stablemate result for genes that are predictive either in the core or in the periphery
load('./case3/r_object/dpt_gene/combine.RData')
mod_combine <- mod

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
mod_core <- change_significance(mod_core)
mod_periphery <- change_significance(mod_periphery)
mod_combine <- change_significance(mod_combine)


############### Figure S8C ############### 
# Make plot
p1 <- plot.SRST2E(mod_core) +  theme(text = element_text(size = 15))+ ylab('') + xlab('Predictivity score: predictive in the tumour core')
p2 <- plot.SRST2E(mod_periphery) +  theme(text = element_text(size = 15)) + ylab('Stability score') + xlab('Predictivity score: predictive in the tumour periphery')
p3 <- plot.SRST2E(mod_combine) +  theme(text = element_text(size = 15)) + ylab('') + xlab('Predictivity score: predictive in the full data')
ggarrange(p2,p3,p1,common.legend = T, ncol = 3, legend = 'right')


############### Figure S9A ############### 
# Stable genes are obtained from the combined analysis
sb <- print(mod_combine)$SB_selected
# Core-specific genes are predictive in the core but unstable
nsb_core <- print(mod_core)$NSB_selected
# Periphery-specific genes are predictive in the periphery but unstable
nsb_periphery<- print(mod_periphery)$NSB_selected
# Stable genes that are found predictive in either locations but unstable are considered false discovery
sb <- sb[!sb %in% nsb_core & !sb %in% nsb_periphery]

# Create column annotation
ano_col <- data.frame(Location = location[order(DPT366)], DPT = DPT366[order(DPT366)])
rownames(ano_col) <- colnames(impdata)[order(DPT366)]

# Heatmap of stable genes
sb_labels <- sb # Only label the genes that are mentioned in the paper
sb_labels[! sb %in% 
            c('CCL3','CCL4','FN1', 'VIM', 'TGFBI', 'LDHA', 'CD83', 'EGR2')] <- ''
smoothed_sb <- apply(impdata[sb,], 1, function(x) predict(gam(x~s(Y, bs = 'cs')))) # Smooth for better visualization
smoothed_sb <- apply(smoothed_sb, 2, function(x) (x-min(x))/(max(x)-min(x)))
pheatmap::pheatmap(t(smoothed_sb[order(Y),]),
                   color = colorRampPalette((brewer.pal(n = 7, name =
                                                          "Blues")))(100),
                   cluster_rows = T, cluster_cols = F, show_colnames = F, 
                   treeheight_row = 20, annotation_col = ano_col,
                   annotation_colors = list(Location = c(Core = 'lavenderblush4',
                                                         Periphery = 'black')), fontsize = 14)

# Heatmap of periphery-specific genes
nsb_periphery_labels <- nsb_periphery # Only label the genes that are mentioned in the paper
nsb_periphery_labels[! nsb_periphery %in% 
                       c('TNF','IL1B','CCL2', 'CSF1')] <- ''
smoothed_nsb_periphery <- apply(impdata[nsb_periphery,], 1, function(x) predict(gam(x~s(Y, bs = 'cs')))) # Smooth for better visualization
smoothed_nsb_periphery <- apply(smoothed_nsb_periphery, 2, function(x) (x-min(x))/(max(x)-min(x)))
pheatmap::pheatmap(t(smoothed_nsb_periphery[order(Y),]),
                   color = colorRampPalette((brewer.pal(n = 7, name =
                                                          "Blues")))(100),
                   cluster_rows = T, cluster_cols = F, show_colnames = F,
                   treeheight_row = 20, annotation_col = ano_col,
                   annotation_colors = list(Location = c(Core = 'lavenderblush4',
                                                         Periphery = 'black')), fontsize = 14)

# Heatmap of core-specific genes
nsb_core_labels <-  nsb_core # Only label the genes that are mentioned in the paper
nsb_core_labels[! nsb_core %in% 
                  c('VCAN','MARCO','HLA-DOA')] <- ''
smoothed_nsb_core <- apply(impdata[nsb_core,], 1, function(x) predict(gam(x~s(Y, bs = 'cs'))))
smoothed_nsb_core <- apply(smoothed_nsb_core, 2, function(x) (x-min(x))/(max(x)-min(x)))
pheatmap::pheatmap(t(smoothed_nsb_core[order(Y),]),
                   color = colorRampPalette((brewer.pal(n = 7, name =
                                                          "Blues")))(100),
                   cluster_rows = T, cluster_cols = F, show_colnames = F,
                   treeheight_row = 20, annotation_col = ano_col,
                   annotation_colors = list(Location = c(Core = 'lavenderblush4',
                                                         Periphery = 'black')), fontsize = 14)


############### Figure S9B (Running S9A is required) ############### 
# Cluster the stable genes
sb_cluster <- hclust(dist(t(smoothed_sb))) # Use the same hierarchical clustering as which used for clustering the rows of the heatmap
sb_cluster <- cutree(sb_cluster, k = 2) # Identify the two major clusters

# Cluster the core-specific genes
nsb_core_cluster <- hclust(dist(t(smoothed_nsb_core))) # Use the same hierarchical clustering as which used for clustering the rows of the heatmap
nsb_core_cluster <- cutree(nsb_core_cluster, k = 2)  # Identify the two major clusters

# Cluster the core-specific genes
nsb_periphery_cluster <- hclust(dist(t(smoothed_nsb_periphery)))
nsb_periphery_cluster <- cutree(nsb_periphery_cluster, k = 2)

# Check cluster id
sb_cluster
# CCL3 and CCL4 negatively correlated with diffusion pseudo time. Therefore, the cluster 1 to which 
# CCL3 and CCL4 belongs is the cluster of negatively stable genes we mentioned in the main text. The cluster2,
# in contrary, is the cluster of positively stable genes.

nsb_core_cluster
# The cluster 1 is the major cluster, thus represents the core-specific genes 

nsb_periphery_cluster
# The cluster 1 is the major cluster, thus represents the periphery-specific genes 

# Now aggregate the gene expression of these major gene modules

# Function for aggregating gene expression, refer to Section 4.3 for details
npca <- function(x){
  pca <- prcomp(x, center = F, scale = F)
  l <- pca$rotation[,1]
  l <- sign(l[which.max(abs(l))]) * l
  l <- l - min(l)
  l <- l/sqrt(sum(l^2))
  as.numeric(x %*% l)
}

# Aggregate expression
neg_stab <- t(impdata)[,sb[sb_cluster == 1]] %>% npca() # Aggregate the negatively stable genes
pos_stab <- t(impdata)[,sb[sb_cluster == 2]] %>% npca() # Aggregate the positively stable genes
core_specific <- t(impdata)[,nsb_core[nsb_core_cluster == 1]] %>% npca() # Aggregate the positively stable genes
periphery_specific <- t(impdata)[,nsb_periphery[nsb_periphery_cluster == 1]] %>% npca() # Aggregate the positively stable genes

# Prepare the data for plotting 
plot_df <- data.frame(expression = c(neg_stab,pos_stab,core_specific, periphery_specific),
                      type = c(rep('Negatively stable',N), rep('Positively stable',N),
                               rep('Core specific',N), rep('Periphery specific',N)),
                      time = c(rep(DPT366, 4)),
                      location = c(rep(location, 4)),
                      simple_location = c(rep(location, 4)))
plot_df$type <- factor(plot_df$type, levels = c('Negatively stable','Positively stable', 
                                                'Periphery specific', 'Core specific')) # Order panels
plot_df$simple_location[plot_df$simple_location %in% c('Tip1','Tip2')] <- 'Core'

# Make the plot
ggplot(data =plot_df ) + geom_point(aes(time,expression, col = location)) + 
  facet_wrap(~type, scales = 'free', ncol = 2) + xlim(c(0,max(DPT366))) + 
  ylim(c(0,NA))+
  scale_color_manual('Location',values = c('lavenderblush4','black')) + 
  ylab('Aggregated expression of gene modules') +
  xlab('DPT') + theme(text = element_text(size = 18), legend.position = 'bottom')  + 
  geom_smooth(aes(time,expression, col = simple_location), method = 'lm')


############### Figure 4B ############### 
p1 <- ggplot() + geom_point(aes(Y, impdata['CCL3',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed CCL3') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['CCL3',], col = location), method = 'lm') + xlab('') + ylim(c(0,NA))
p2<- ggplot() + geom_point(aes(Y, impdata['CCL4',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed CCL4') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['CCL4',], col = location), method = 'lm') + xlab('DPT')+ ylim(c(0,NA))
ggarrange(p1,p2, common.legend = T, ncol = 1, legend = 'none', align = 'hv')


############### Figure 4C ############### 
p1 <- ggplot() + geom_point(aes(Y, impdata['TNF',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed TNF') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['TNF',], col = location), method = 'lm') + xlab('') + ylim(c(0,NA))
p2<- ggplot() + geom_point(aes(Y, impdata['CCL2',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed CCL2') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['CCL2',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p3<- ggplot() + geom_point(aes(Y, impdata['IL1B',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed IL1B') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['IL1B',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p4<- ggplot() + geom_point(aes(Y, impdata['CSF1',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed CSF1') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['CSF1',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
ggarrange(p1,p2,p3,p4, common.legend = T, ncol = 2, legend = 'none', align = 'hv')


############### Figure 4E ############### 
p1 <- ggplot() + geom_point(aes(Y, impdata['MARCO',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed MARCO') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['MARCO',], col = location), method = 'lm') + xlab('') + ylim(c(0,NA))
p2<- ggplot() + geom_point(aes(Y, impdata['HLA-DOA',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed HLA-DOA') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['HLA-DOA',], col = location), method = 'lm') + xlab('DPT')+ ylim(c(0,NA))
ggarrange(p1,p2, common.legend = T, ncol = 1, legend = 'none', align = 'hv')


############### Figure 4F ############### 
ggplot() + geom_point(aes(Y, impdata['VCAN',],col = location))  + 
  scale_color_manual('Location', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed VCAN') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['VCAN',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))+
  xlab('DPT')


############### Figure S8D ############### 
p1 <- ggplot() + geom_point(aes(Y, impdata['EGR2',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed EGR2') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['EGR2',], col = location), method = 'lm') + xlab('DPT') + ylim(c(0,NA))
p2<- ggplot() + geom_point(aes(impdata['CCL3',], impdata['EGR2',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed EGR2') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(impdata['CCL3',], impdata['EGR2',], col = location), method = 'lm') + xlab('Sincast imputed CCL3')+ ylim(c(0,NA))
ggarrange(p1,p2, common.legend = T, ncol = 1, legend = 'none', align = 'hv')


############### Figure S8E ############### 
p1 <- ggplot() + geom_point(aes(Y, impdata['HLA-DRA',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed HLA-DRA') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['HLA-DRA',], col = location), method = 'lm') + xlab('') + ylim(c(0,NA))
p2<- ggplot() + geom_point(aes(Y, impdata['HLA-DOA',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed HLA-DOA') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['HLA-DOA',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p3<- ggplot() + geom_point(aes(Y, impdata['HLA-DMA',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed HLA-DMA') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['HLA-DMA',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p4<- ggplot() + geom_point(aes(Y, impdata['HLA-DPA1',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed HLA-DPA1') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['HLA-DPA1',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p5<- ggplot() + geom_point(aes(Y, impdata['HLA-DQA1',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed HLA-DQA1') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['HLA-DQA1',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p6<- ggplot() + geom_point(aes(Y, impdata['HLA-DQA2',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed HLA-DQA2') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['HLA-DQA2',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))

ggarrange(p1,p2,p3,p4,p5,p6, common.legend = T, ncol = 3, legend = 'none', align = 'hv')


############### Figure S8F ############### 
p1 <- ggplot() + geom_point(aes(Y, impdata['TMEM119',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed TMEM119') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['TMEM119',], col = location), method = 'lm') + xlab('') + ylim(c(0,NA))
p2<- ggplot() + geom_point(aes(Y, impdata['P2RY12',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed P2RY12') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['P2RY12',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p3<- ggplot() + geom_point(aes(Y, impdata['GPR34',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed GPR34') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['GPR34',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
p4<- ggplot() + geom_point(aes(Y, impdata['OLFML3',],col = location))  + 
  scale_color_manual('Sample type', values = c('lavenderblush4','black')) +
  ylab('Sincast imputed OLFML3') + theme(text = element_text(size= 15)) + 
  geom_smooth(aes(Y, impdata['OLFML3',], col = location), method = 'lm') + xlab('')+ ylim(c(0,NA))
ggarrange(p1,p2,p3,p4, common.legend = T, ncol = 2, legend = 'none', align = 'hv')


############### FOR THE FOLLOWING PLOT, CLEAN UP MEMORY, RESTART R AND RERUN PRELOADING ############### 


############### The following preloading was used for Figure 4D (MUST RUN) ############### 
require(destiny)
require(dplyr)
require(qgraph)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate/')
source('./stablemate.R')

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


############### Figure 4D ############### 
# Get stable genes
# Stable genes that are found predictive in either locations but unstable are considered false discovery

# Get stable genes for predicting CCL3
load('./case3/r_object/ccl3_gene/combine.RData')
mod <- change_significance(mod)
sb_ccl3 <- print(mod)$SB_selected
load('./case3/r_object/ccl3_gene/core.RData')
mod <- change_significance(mod)
sb_ccl3 <- sb_ccl3[!sb_ccl3 %in% print(mod)$NSB_selected] 
load('./case3/r_object/ccl3_gene/periphery.RData')
mod <- change_significance(mod)
sb_ccl3 <- sb_ccl3[!sb_ccl3 %in% print(mod)$NSB_selected]

# Get stable genes for predicting CCL4
load('./case3/r_object/ccl4_gene/combine.RData')
mod <- change_significance(mod)
sb_ccl4 <- print(mod)$SB_selected
load('./case3/r_object/ccl4_gene/core.RData')
mod <- change_significance(mod)
sb_ccl4 <- sb_ccl4[!sb_ccl4 %in% print(mod)$NSB_selected]
load('./case3/r_object/ccl4_gene/periphery.RData')
mod <- change_significance(mod)
sb_ccl4 <- sb_ccl4[!sb_ccl4 %in% print(mod)$NSB_selected]

# Get stable genes for predicting TNF
load('./case3/r_object/tnf_gene/combine.RData')
mod <- change_significance(mod)
sb_tnf <- print(mod)$SB_selected
load('./case3/r_object/tnf_gene/core.RData')
mod <- change_significance(mod)
sb_tnf <- sb_tnf[!sb_tnf %in% print(mod)$NSB_selected]
load('./case3/r_object/tnf_gene/periphery.RData')
mod <- change_significance(mod)
sb_tnf <- sb_tnf[!sb_tnf %in% print(mod)$NSB_selected]

# Get stable genes for predicting CCL2
load('./case3/r_object/ccl2_gene/combine.RData')
mod <- change_significance(mod)
sb_ccl2 <- print(mod)$SB_selected
load('./case3/r_object/ccl2_gene/core.RData')
mod <- change_significance(mod)
sb_ccl2 <- sb_ccl2[!sb_ccl2 %in% print(mod)$NSB_selected]
load('./case3/r_object/ccl2_gene/periphery.RData')
mod <- change_significance(mod)
sb_ccl2 <- sb_ccl2[!sb_ccl2 %in% print(mod)$NSB_selected]

# Get stable genes for predicting IL1B
load('./case3/r_object/il1b_gene/combine.RData')
mod <- change_significance(mod)
sb_il1b <- print(mod)$SB_selected
load('./case3/r_object/il1b_gene/core.RData')
mod <- change_significance(mod)
sb_il1b <- sb_il1b[!sb_il1b %in% print(mod)$NSB_selected]
load('./case3/r_object/il1b_gene/periphery.RData')
mod <- change_significance(mod)
sb_il1b <- sb_il1b[!sb_il1b %in% print(mod)$NSB_selected]

# Get stable genes for predicting CSF1
load('./case3/r_object/csf1_gene/combine.RData')
mod <- change_significance(mod)
sb_csf1 <- print(mod)$SB_selected
load('./case3/r_object/csf1_gene/core.RData')
mod <- change_significance(mod)
sb_csf1 <- sb_csf1[!sb_csf1 %in% print(mod)$NSB_selected]
load('./case3/r_object/csf1_gene/periphery.RData')
mod <- change_significance(mod)
sb_csf1 <- sb_csf1[!sb_csf1 %in% print(mod)$NSB_selected]

# Get stable genes for predicting diffusion pseudo time 
# The diffusion pseudo time starting from the 366th cell (the cell with the smallest 
# DC1 value) is chosen as the response for stablemate analysis.
DPT366 <- DPT$DPT366
load('./case3/r_object/dpt_gene/combine.RData')
mod <- change_significance(mod)
sb_dpt366 <- print(mod)$SB_selected
load('./case3/r_object/dpt_gene/core.RData')
mod <- change_significance(mod)
sb_dpt366 <- sb_dpt366[!sb_dpt366 %in% print(mod)$NSB_selected]
load('./case3/r_object/dpt_gene/periphery.RData')
mod <- change_significance(mod)
sb_dpt366 <- sb_dpt366[!sb_dpt366 %in% print(mod)$NSB_selected]

# Get genes that are used as nodes of the network 
nodes <- c('DPT','CCL3','CCL4','TNF','CCL2','IL1B','CSF1') # The major node
nodes <- c(nodes, sb_ccl3, sb_ccl4, sb_tnf, sb_ccl2,
           sb_il1b, sb_csf1, sb_dpt366) # Stable genes of each run of stablemate are also included as minor nodes
nodes <- unique(nodes) # Remove duplication
n_nodes <- length(nodes)

# Load KEGG cytokines, used for annotate nodes

# Load KEGG interleukins
interleukins <- readLines('./case3/downloaded_file/kegg_cytokines/interleukins')
# Only keep interleukins that are in the query data, same for the rest of cytokines
interleukins <- intersect(rownames(impdata),interleukins) 

# Load KEGG interferons
interferons <- readLines('./case3/downloaded_file/kegg_cytokines/interferons')
interferons <- intersect(rownames(impdata),interferons) 

# Load KEGG colony Stimulating factors
csfs <- readLines('./case3/downloaded_file/kegg_cytokines/csfs')
csfs <- intersect(csfs,rownames(impdata))

# Load KEGG chemokines
chemokines <-  readLines('./case3/downloaded_file/kegg_cytokines/chemokines')
chemokines <- intersect(chemokines,rownames(impdata))

# Load KEGG tumor necrosis factor
tnfs <-  readLines('./case3/downloaded_file/kegg_cytokines/tnfs')
tnfs <- intersect(tnfs,rownames(impdata))

# Load KEGG tumor growth factor
tgfs <-  readLines('./case3/downloaded_file/kegg_cytokines/tgfs')
tgfs <- intersect(tgfs,rownames(impdata))

# Load KEGG growth factors
gfs <-  readLines('./case3/downloaded_file/kegg_cytokines/gfs')
gfs <- intersect(gfs,rownames(impdata))

# Combine cytokines
cytokines <- list(interleukins = interleukins, interferons= interferons, csfs =csfs, 
                  chemokines = chemokines, tnfs = tnfs, tgfs = tgfs, gfs = gfs)

# Make the network

# Create the adjacency matrix of the network
adj <- matrix(0, nrow = n_nodes, ncol = n_nodes, 
              dimnames =  list(nodes, nodes))
# Genes that are stable predictors of CCL3 are connected with CCL3. Same for the other major nodes
adj['CCL3',sb_ccl3] <- 1 
adj['CCL4',sb_ccl4] <- 1
adj['TNF',sb_tnf] <- 1
adj['CCL2',sb_ccl2] <- 1
adj['IL1B',sb_il1b] <- 1
adj['CSF1',sb_csf1] <- 1
adj['DPT',sb_dpt366] <- 1
adj <- (adj | t(adj))

# Calculate correlation between the nodes
cor_mat <- as.matrix(rbind(DPT366, impdata[nodes[-1],])) # Concatenate diffusion pseudi time and gene expression
cor_mat <- cor(t(cor_mat)) # Calculate correlation

# Annotate nodes by their cytokine type
node_type <- c()

convert_zero_length_to_others <- function(x) ifelse(length(x) == 0, 'others', x)  

for(i in 1:n_nodes){
  
  # Check each node belongs to which cytokine category
  node_type[i] <- sapply(cytokines, function(c) nodes[i] %in% c) %>% which() %>% names() %>% 
    convert_zero_length_to_others()
  
}

# Highlight non-cytokine nodes that are selected as stable predictors of at least two major nodes 
# ('CCL3','CCL4','TNF','CCL2','IL1B','CSF1')
node_type[node_type == 'others' &
            rowSums(adj) >= 2 &
            1:ncol(adj) > 1] <- 'noted'

# Annotate node labels by their cytokine type
label_type <- node_type
label_type[1:7] <- 'noted' # To note all major nodes

# Create color schemes for labels and nodes
node_cols <- c('noted' = 'black','others' = 'lavenderblush4', 'chemokines' = 'red', 'csfs' = 'blue', 'interleukins' = 'darkgreen', 'interferons' = 'orangered3',
               'tnfs' = 'purple', 'tgfs' = 'darkolivegreen4', 'gfs' = 'deeppink4')
label_cols <- c('noted' = 'black','others' = 'lavenderblush4', 'chemokines' = 'red', 'csfs' = 'blue', 'interleukins' = 'darkgreen', 'interferons' = 'orangered3',
                'tnfs' = 'purple', 'tgfs' = 'darkolivegreen4', 'gfs' = 'deeppink4')

# Generate color vectors for labels and nodes
node_col <- alpha(node_cols[node_type], alpha = 0.8)
label_col <- alpha(label_cols[label_type], alpha = 0.8)

# Generate size vectors for labels and nodes
node_size <-  rep(0.5, n_nodes)
node_size[1:7] <- 6 # The major nodes are larger, whereas the minor nodes are merely invisible

label_size <- rep(1, n_nodes)
label_size[rowSums(adj) >1] <- 13
label_size[label_type != 'others'] <- 13
label_size[1:7] <- 1

# Wight adjacency matrix by correlation
adj_cor <- adj * cor_mat

# Create spring layout of the network based on unweighted adjacency
g <- qgraph(adj, groups = node_type, vsize = node_size, label.cex = label_size, label.color = label_col,
            color = node_col, esize=2,  repulsion = 1, labels = rownames(adj), borders = FALSE, edge.color = 'grey',
            layout = 'spring')

# Based on the created layout above, regenerate the network but with edges weighted and colored by correlation 
qgraph(adj_cor, groups = node_type, vsize = node_size, label.cex = label_size, label.color = label_col,
       color = node_col, esize=2,  repulsion = 1, labels = rownames(adj), borders = FALSE, 
       layout = g$layout)
