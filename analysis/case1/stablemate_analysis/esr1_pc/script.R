# Load stablemate functions
# Some essential package has already been loaded
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate')
source('./stablemate.R')

############### Parallel computing setting, to be ran in slurm system ############### 
cl <- startMPIcluster(verbose=TRUE) # By default will start one fewer slave
# Than elements in tasks requested (to account for master)
registerDoMPI(cl)
niter <- 20


############### Data preprocessing ############### 
load('./case1/data/tcga/TCGA-BRCA.RData')

#Filter for protein coding genes
data <- data[data@rowRanges$gene_type == 'protein_coding',]
data <- data[!duplicated(data@rowRanges$gene_name),]
rownames(data) <- data@rowRanges$gene_name

#Set response, predictors, and environments
gene <- 'ESR1'
X_full <- log(t(data@assays@data$tpm_unstrand)+1) # log-tmp transformation
rownames(X_full) <- colnames(data)
colnames(X_full) <- rownames(data)

Y_full <- X_full[,gene] # Remove ESR1
X_full <- X_full[,colnames(X_full)!=gene]  
meta_data_full <- data@colData

#Get ER+ and normal samples
idx <- meta_data_full$paper_BRCA_Subtype_PAM50 %in% c('LumA','LumB') | 
  meta_data_full$sample_type == 'Solid Tissue Normal'
X <- X_full[idx ,]
meta_data <- meta_data_full[idx,]
Y <- Y_full[idx]
envs <- meta_data$sample_type

#Free memory
gene_anno <- data@rowRanges
rm(data)


############### PCA ############### 
pca <- prcomp(X)

# Find the elbow point for choosing the number of PCs to use
find.elbow <- function(x, y){
  
  n <- length(x)
  scaled.y <- y * max(x)/max(y)
  r <- -(scaled.y[1] - scaled.y[n])/(x[1]-x[n])
  c <- - scaled.y[1] - r * x[1]
  best.index <- which.max(abs(r * x +  scaled.y + c))
  best.index
  
}

elbow_point <- find.elbow(1:length(pca$sdev^2),pca$sdev^2)

############### Run stablemate ############### 
mod_list <- foreach(i=1:niter, .verbose = T) %dopar% {
  library(progress)
  library(arm)
  library(Matrix)
  SRST2E(Y, pca$x[,1:elbow_point], envs, K = 100,
         return_prob = F,subsample_method = 'all')
}
mod <- combine_mc(mod_list)
save.image('./case1/r_object/esr1_pc/ESR1.RData')