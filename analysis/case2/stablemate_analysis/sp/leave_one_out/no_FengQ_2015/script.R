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
# Prepare the input for stablemate
data_list <- readRDS('./case2/data/processed_data_sp.RDS')
envs <- data_list[[3]]$cohort
Y <- as.numeric(factor(data_list[[3]]$sample_type))-1
X <- data_list[[1]]

# Split training and testing
Ytrain <- Y[envs != "FengQ_2015"]
Xtrain <- X[envs != "FengQ_2015" , ]
Ytest <- Y[envs == "FengQ_2015"]
Xtest <- X[envs == "FengQ_2015" , ]
envs <- envs[envs != "FengQ_2015"]


############### Run stablemate ############### 
K <- 100
mod_list <- foreach(i=1:niter, .verbose = T) %dopar% {
  library(progress)
  library(arm)
  library(Matrix)
  mod <- SRST2E(Ytrain, Xtrain, envs, K = K, b = 100/ncol(Xtrain),         
                pred_st2_control = list(obj_fun = obj_BIC_logit, reg_fun = reg_logit),
                stab_st2_control = list(obj_fun = obj_PSS_logit, reg_fun = reg_logit),
                return_prob = F,subsample_method = 'all')
}
mod <- combine_mc(mod_list)
save.image('./case2/r_object/sp/leave_one_out/no_FengQ_2015.RData')
