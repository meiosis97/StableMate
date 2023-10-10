# Load stablemate functions
# Some essential package has already been loaded
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate')
source('./stablemate.R')
require(destiny)

############### Parallel computing setting, to be ran in slurm system ###############
cl <- startMPIcluster(verbose=TRUE) # By default will start one fewer slave
# Than elements in tasks requested (to account for master)
registerDoMPI(cl)
niter <- 20


############### Data preprocessing ###############
# Load the Seurat object storing the raw, the log-transformed data and the metadata
seurat_obj <- readRDS('./case3/data/processed_query_data/seurat_data.rds')  
impdata <- readRDS('./case3/data/processed_query_data/impdata.rds') # Load the Sincast imputed data
metadata <- seurat_obj@meta.data # Extract the metadata from the Seurat object

# Set response, predictors, and environments
envs <- metadata$Cluster_2d
envs[envs ==7] <- 'Periphery'
envs[envs =='8'] <- 'Core'

Y <- impdata['IL1B',]
X <- as.matrix(t(impdata[rownames(impdata) != 'IL1B',]))
Y_sub <- Y[envs == 'Periphery']
X_sub <- X[envs == 'Periphery',]

# Define a special objective function, so that the predictivity of predictors are measured within the periphery
obj_BIC <- function(Y, X, env, var_in, sample_id){
  
  sample_id = which(env == 'Periphery') # Only consider samples within the periphery
  n_sample <- length(sample_id)
  
  # Weight samples by environment sizes
  w <- 1/(table(env[sample_id])[env[sample_id]])
  w <- w/sum(w)*n_sample
  n_in <- sum(var_in)
  if(n_in>0){
    mod <- lm(Y[sample_id]~., data.frame(X[sample_id, var_in==1, drop = FALSE]), weights = w)
  }else{
    mod <- lm(Y[sample_id]~1, weights = w)
  }
  (mod$rank+1) * log(n_sample) -
    2*sum(w*dnorm(Y[sample_id], fitted(mod), sqrt(sum(residuals(mod)^2)/n_sample), log =T))
  
}


############### Prefiltering of predictor sets using stability selection (see Supplementary Methods for detail description) ###############
K <- 100 # Ensemble size
size <- 100 # Prefilter size 
N <- length(Y_sub)

# Lasso selection with parallel computing
lasso_list <- foreach(i=1:niter, .verbose = T) %dopar% {
  
  require(glmnet)
  lasso_sel <- matrix(0, ncol = ncol(X_sub), nrow = K)
  colnames(lasso_sel) <- colnames(X_sub)
  
  for(j in 1:K){
    
    idx <- sample(1:N, ceiling(N/2)) # Sample data in half
    n_sample <- length(idx)
    w <- 1/(table(envs[idx])[envs[idx]]) # Weight samples by sizes of their environments
    w <- w/sum(w)*n_sample
    mod <- glmnet(X_sub[idx,],Y_sub[idx], weights = w) # Run Lasso on subsampled data
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
save.image('./case3/r_object/il1b_gene/periphery.RData')