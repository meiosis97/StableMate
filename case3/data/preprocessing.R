############### The following preloading was used for the entire script (MUST RUN) ############### 
############### All the plot sessions can be ran independently of each others after running the preloading ############### 
require(Seurat)
require(ggplot2)
require(Sincast)
require(destiny)
require(plotly)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate/')


############### Select the top 2000 most variable genes using Seurat ############### 
# Load data
data <- read.table('./case3/data/query_data/GBM_raw_gene_counts.csv')
metadata <- read.table('./case3/data/query_data/GBM_metadata.csv')

# Clean up cell names
colnames(data) <- substr(colnames(data), start = 2, stop = nchar(colnames(data)))
# Only retain myeloid cells, which correspond to cluster 7 and 8 of the original annotation
data <- data[, metadata$Cluster_2d %in% c(7,8)]
metadata <- subset(metadata, metadata$Cluster_2d %in% c(7,8))

# Create a Seurat object without quality control
seurat_obj <- CreateSeuratObject(counts = data, meta.data = metadata,min.cells = 0, min.features = 0)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000) # Log-normalization
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000) # Find the most variable genes

# Save the data as a Seurat object
saveRDS(seurat_obj, file = './case3/data/processed_query_data/seurat_data.rds')


############### Data visualization (not necessary, require the previous code section to be ran first) ############### 
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)) # Data scaling (centering)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj)) # PCA
seurat_obj <- RunTSNE(seurat_obj, dims = 1:50) # TSNE
seurat_obj <- RunUMAP(seurat_obj, dims = 1:50) # UMAP
DimPlot(seurat_obj, reduction = "pca", group.by = 'Cluster_2d')
DimPlot(seurat_obj, reduction = "tsne", group.by = 'Cluster_2d')
DimPlot(seurat_obj, reduction = "umap", group.by = 'Cluster_2d')


############################ Sincast imputation ################################
# Load the Seurat object
seurat_obj <- readRDS('./case3/data/processed_query_data/seurat_data.rds')

# Store the gene expression data of the variable genes as an sce object and run Sincast imputation
vf <- VariableFeatures(seurat_obj) # Extract the variable genes
sce_vf <- createSce(data = seurat_obj@assays$RNA@data[vf,], colData = seurat_obj@meta.data) # Create the sce object
sce_vf <- sincastImp(sce_vf, dologScale = F, col.by = 'Cluster_2d') # Sincast imputation

# Based on the imputation operator learned on the variable genes, impute the full data.
impdata_full <- randomWalk(sce_vf, seurat_obj@assays$RNA@data, logScale = F, t= 3)

# Create an sce object to store the full data before and after imputation
sce_full <- SingleCellExperiment(list(before = seurat_obj@assays$RNA@data,
                                      after =  impdata_full))

# Post imputation scaling to avoid overfitting
sce_full <- postScale(query = sce_full, preImpAssay = 'before',
                      postImpAssay = 'after', dologScale = F)

# Save the imputed full data 
saveRDS(sce_full@assays@data$SincastScaledData,  './case3/data/processed_query_data/impdata.rds')


############################ Diffusion map and diffusion pseudo time ################################
# Load the Seurat object
seurat_obj <- readRDS('./case3/data/processed_query_data/seurat_data.rds')

# Run diffusion map on the log-normalized gene expression data of the variable genes
vf <- VariableFeatures(seurat_obj)
DM <- DiffusionMap(data = t(as.matrix(seurat_obj@assays$RNA@data[vf,])), n_pcs = 50)

# Diffusion pseudo time learning
DPT <- DPT(DM)

# Visualize diffusion pseudo time
# Rename cluster id to cell locations, which will be used to annotate query cells
location <- seurat_obj@meta.data$Cluster_2d 
location[location == 7] <- 'Peripheral'
location[location =='8'] <- 'Core'

# 2D ggplot
ggplot() + geom_point(aes(DM$DC1, DM$DC2,
                          col = location)) +
  xlab('DC1') + ylab('DC2') + 
  scale_color_manual('Location', values = c('lavenderblush4','black')) +
  theme(text = element_text(size = 15))

# 3D plotly
plot_ly(data = data.frame(attr(DPT@dm,'eigenvectors'))) %>% 
  add_markers(x=~DC1, y =~DC2, z = ~DC3, 
              color = location,
              colors = c('lavenderblush4','black'))

# A clear cell state transition is observed. Consider the cell with the smallest DC1 coordination as 
# the root of transition
which.min(DM$DC1)
# Cell 366 is set as the root.

# Rerun diffusion pseudo time learning, but with further specification on the tips of cell trajectory
DPT <- DPT(DM, tips = c(which.min(DM$DC1), which.max(DM$DC3),which.min(DM$DC3)))

# Save the diffusion pseudo time object
saveRDS(DPT, file = './case3/data/processed_query_data/dpt.rds')
