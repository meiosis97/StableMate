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
data_list <- readRDS('./case2/data/processed_data_pw.RDS')

# Make environment variable. Samples outside of the target cohort are labelled 'others'
envs <- data_list[[3]]$cohort
envs[envs != "study"] <- 'others'

Y <- as.numeric(factor(data_list[[3]]$sample_type))-1
X <- data_list[[1]]

# Define a special objective function, so that the predictivity of predictors are measured within the study cohort.
obj_BIC_logit <- function(Y, X, env, var_in, sample_id = 1:length(Y)){
  
  sample_id <- which(env == "study") # Only consider samples within the study cohort.
  n_sample <- length(sample_id) 
  
  # Weight samples by cohort-case (disease status) sizes
  agg_group <- paste(Y, env)
  w <- 1/(table(agg_group[sample_id])[agg_group[sample_id]])
  w <- w/sum(w)*n_sample
  
  n_in <- sum(var_in)
  if(n_in>0){
    mod <- bayesglm(Y[sample_id]~., data.frame(X[sample_id, var_in==1, drop = FALSE]),
                    family = 'binomial', weights = w, control = glm.control(maxit = 200))
  }else{
    mod <- bayesglm(Y[sample_id]~1, family = 'binomial',
                    weights = w, control = glm.control(maxit = 200))
  }
  (mod$rank) * log(n_sample) - 2 * sum(w*dbinom(Y[sample_id], 1, fitted(mod), log =T))
  
}


############### Run stablemate ############### 
K <- 100
mod_list <- foreach(i=1:niter, .verbose = T) %dopar% {
  library(progress)
  library(arm)
  library(Matrix)
  mod <- SRST2E(Y, X, envs, K = K, max_size = sum(envs == "study"),
                pred_st2_control = list(obj_fun = obj_BIC_logit, reg_fun = reg_logit),
                stab_st2_control = list(obj_fun = obj_PSS_logit, reg_fun = reg_logit),
                return_prob = F,subsample_method = 'all')
}
mod <- combine_mc(mod_list)
save.image('./case2/r_object/pw/one_by_one/study.RData')
