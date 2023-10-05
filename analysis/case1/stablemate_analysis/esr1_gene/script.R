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
  
# Set response, predictors, and environments
gene <- 'ESR1'
X_full <- log(t(data@assays@data$tpm_unstrand)+1) # log-tmp transformation
rownames(X_full) <- colnames(data)
colnames(X_full) <- rownames(data)

Y_full <- X_full[,gene] # Remove ESR1
X_full <- X_full[,colnames(X_full)!=gene]  
meta_data_full <- data@colData

# Get ER+ and normal samples
idx <- meta_data_full$paper_BRCA_Subtype_PAM50 %in% c('LumA','LumB') | 
        meta_data_full$sample_type == 'Solid Tissue Normal'
X <- X_full[idx ,]
meta_data <- meta_data_full[idx,]
Y <- Y_full[idx]
envs <- meta_data$sample_type

# Free memory
gene_anno <- data@rowRanges
rm(data)

############### Prefiltering of predictor sets using stability selection (see Supplementary Methods for detail description) ###############
K <- 100 # Ensemble size
size <- 100 # Prefilter size 
N <- length(Y)

# Lasso selection with parallel computing
lasso_list <- foreach(i=1:niter, .verbose = T) %dopar% {
  
  require(glmnet)
  lasso_sel <- matrix(0, ncol = ncol(X), nrow = K)
  colnames(lasso_sel) <- colnames(X)
  
  for(j in 1:K){
    
    idx <- sample(1:N, ceiling(N/2)) # Sample data in half
    n_sample <- length(idx)
    w <- 1/(table(envs[idx])[envs[idx]]) # Weight samples by sizes of their environments
    w <- w/sum(w)*n_sample
    mod <- glmnet(X[idx,],Y[idx], weights = w) # Run Lasso on subsampled data
    sel.matrix <- (mod$beta != 0) # Indicator matrix that record lasso selection path
    first.entrance <- apply(sel.matrix, 1, which.max) # Calculate at which lambda each predictor is first selected
    first.entrance[which(apply(sel.matrix, 1, sum) == 0)] <- Inf # If not selected in the entire path, set the selection order to infinity
    num_noinf <- sum(first.entrance != Inf) # Calculate the number of predictors that are not selected by the entire selection path
    selected_var <- names(which(first.entrance <= 
                                  sort(first.entrance)[min(size,num_noinf)])) # Get the name of first 100 predictors entering the path
    lasso_sel[j,selected_var] <- 1
    
  }
  
  Matrix(lasso_sel)
  
}

# Prepare the predictor pool for stablemate
lasso_sel <- do.call(rbind,lasso_list)
selected_var <- names(which(colSums(lasso_sel)!=0))
rm(lasso_sel)
lasso_list <- lapply(lasso_list, function(x) cbind('pseudo_variable' = 1,
                                                   x[,selected_var, drop = FALSE]))
X <- X[,selected_var]


############### Run stablemate ############### 
mod_list <- foreach(i=1:niter, .verbose = T) %dopar% {
library(progress)
library(arm)
library(Matrix)
  SRST2E(Y, X, envs, K = K, pred_st2_control = list(var_pool = lasso_list[[i]]),
         return_prob = F,subsample_method = 'all')
}
mod <- combine_mc(mod_list)
save.image('./case1/r_object/esr1_gene/ESR1.RData')