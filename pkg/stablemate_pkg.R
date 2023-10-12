############### (This needs to be uncommented if haven't done the installation) ##########
# Packages needs to be installed
# required_pkgs <- c('Matrix', 'dplyr', 'doParallel', 'doMPI',
#                    'ggplot2', 'reshape2', 'ggrepel', 'MASS', 'progress', 'arm')
# to_install <- required_pkgs[!required_pkgs %in% installed.packages()[,1]]
# install.packages(to_install)

# The following functions will be removed in a proper package
`%>%` <- getFromNamespace("%>%", "magrittr")

# Generator of beta-binomial distribution
rbeta_binom <- function(n, size, shape1, shape2) rbinom(n, size, prob = rbeta(n, shape1, shape2))


# Probability mass function of beta-binomial distribution
dbeta_binom <- function(x, size, shape1, shape2) choose(size,x) * beta(x + shape1, size - x + shape2) / beta(shape1, shape2)


# The main ST2E function for stochastic stepwise variable selection
st2e <- function(Y, X, env = NULL, 
                 obj_fun = obj_bic, sel_fun = sel_min, reg_fun = reg_ols,
                 pred_pool = NULL, start_set = 'random_start', max_size = NULL, 
                 alpha = 1, beta = 5, lambda = 0.5, 
                 a = 1, b = 1, t = 0, K = 100, 
                 do_switch = FALSE, 
                 ret_mod = TRUE,
                 ret_imp = TRUE, calc_imp_ctrl = NULL, 
                 ncore = NULL, par_method = c('SNOW','MC','MPI'), verbose = TRUE,
                 chunk_size = ceiling(K/ncore), fun_export = NULL, pkg_export = NULL,
                 drop_pred = TRUE){
  
  ############### Initialize ############### 
  N <- length(Y)
  P <- ncol(X) + 1
  
  # If environment variable is not provided, create a pseudo environment that includes all samples
  if(is.null(env)) env <- rep('pseudo_environment', N)
  env <- as.character(env)
  pred_names <- c('pseudo_predictor', colnames(X))
  scores <- c()  # A vector for storing objective scores
  models <- list() # A list for storing fitted models
  
  # Redefine objective function so that when the the new objective function does not consider the inclusion of pseudo-predictor
  new_obj_fun <- function(Y, X, env, pred_in) obj_fun(Y, X, env, pred_in[-1])
  
  
  ############### Initialize predictor pools ############### 
  # By default, the predictor pool is set to the column names of X plus a pseudo-predictor (all predictors are considered)
  if(is.null(pred_pool)){
    pred_pool <- Matrix::Matrix(1, nrow = K, ncol = P, dimnames = list(NULL, pred_names))
    
    # If a vector of predictor names is provided to the pred_pool argument, convert it to an indicator matrix
  }else if(is.vector(pred_pool)){ 
    # The predictor pool must be a subset of all predictors in X
    if(!all(pred_pool %in% pred_names)) stop('The predictor pool provided is not a subset of the column names of X.')
    tmp <- Matrix::Matrix(0, nrow = K, ncol = P, dimnames = list(NULL, pred_names))
    tmp[,pred_pool] <- 1 
    pred_pool <- tmp
    
    # If an indicator matrix is provided to the pred_pool argument, check whether the format is correct
  }else{
    # If predictor pools provided are the results of StableMate pre-filtering based on Lasso, scale the
    # importance score of the pseudo-predictor
    if(!is.null(attributes(pred_pool)$stbm_lasso) & is.null(calc_imp_ctrl$scale_psd_imp)){
      calc_imp_ctrl$scale_psd_imp <- TRUE
      
    }
    
    # The number of predictor pools provided must equal to the size of the ensemble
    if(K != nrow(pred_pool)){
      message('nrow of pred_pool does not equal to K. Resetting K to nrow of pred_pool')
      K <- nrow(pred_pool)
    }
    # The predictor pools must be subsets of all predictors in X
    if(!all(colnames(pred_pool) %in% pred_names)) stop('The predictor pool provided is not a subset of the column names of X.')
    tmp <- Matrix::Matrix(0, nrow = K, ncol = P, dimnames = list(NULL, pred_names))
    tmp[,colnames(pred_pool)] <- pred_pool
    pred_pool <- tmp
    
  }
  
  # Remove predictors that are not in the predictor pool to save space and time
  if(drop_pred){
    in_pool <- c(TRUE, (colSums(pred_pool) > 0)[-1]) # Always keep the pseudo-predictor
    P <- sum(in_pool)  
    X <- X[,in_pool[-1]]
    pred_pool <- pred_pool[,in_pool]
    pred_names <- colnames(pred_pool)
  }
  
  
  ############### Initialize starting sets based on which ST2 search starts ############### 
  # By default, start from random predictor sets
  do_rand_start <- FALSE
  
  # By default, the starting sets are set as empty
  if(is.null(start_set)){
    start_set <- Matrix::Matrix(0, nrow = K, ncol = P, dimnames = list(NULL, pred_names))
    
    # If a vector is provided to the start_set argument, convert it to an indicator matrix.
  }else if(is.vector(start_set)){ 
    
    # If the starting sets are random, first initialize them as empty sets and perform forward selection with special tuning
    if(length(start_set) == 1 & start_set[1] == 'random_start'){
      start_set <-  Matrix::Matrix(0, nrow = K, ncol = P, dimnames = list(NULL, pred_names))
      do_rand_start <- TRUE
    }else{
      # The starting set must be a subset of all predictors in X
      if(!all(start_set %in% pred_names)) stop('The starting set provided is not a subset of the column names of X.')
      tmp <-  Matrix::Matrix(0, nrow = K, ncol = P, dimnames = list(NULL, pred_names))
      tmp[,start_set] <- 1
      start_set <- tmp
      # The starting set must be a subset of all the predictor pools
      if(any(start_set - pred_pool == 1)) stop('The starting set provided is not a subset of at least one predictor pools')
    }
    
    # If an indicator matrix is provided to the start_set argument, check whether the format is correct.
  }else{
    # The number of starting sets provided must equal to the size of the ensemble
    if(K != nrow(start_set)) stop('nrow of start_set does not equal to nrow of pred_pool')
    # The starting sets must be subsets of all predictors in X
    if(!all(colnames(start_set) %in% pred_names)) stop('The starting set provided is not a subset of the column names of X.')
    tmp <-  Matrix::Matrix(0, nrow = K, ncol = P, dimnames = list(NULL, pred_names))
    tmp[,colnames(start_set)] <- start_set
    start_set <- tmp
    # All starting sets must be the subset of their corresponding predictor pool
    if(any(start_set-pred_pool == 1)) stop('At least one starting set is not the subset of its corresponding predictor pool.')
    
  }
  
  ensemble <- as.matrix(start_set) # Initialize the ensemble 
  if(do_rand_start) start_set <- as.matrix(start_set)
  
  ############### Initialize the maximum allowed model sizes ############### 
  # Flagging whether pseudo-predictors are added into the predictor pools
  
  if(sum(pred_pool[,1])>0){
    message('Pseudo-predictors are added into the predictor pools')
    psd_pred <- TRUE
    
  }else{
    message('Pseudo-predictors are not added into the predictor pools')
    psd_pred <- FALSE
    
  } 
  
  # Calculate the sizes of the predictor pools
  pool_size <- rowSums(pred_pool)
  
  # If the maximum sizes of selections (the max_size argument) are not given, by default set them to the sizes of the predictor pools.
  if(is.null(max_size)){
    max_size <- pool_size
    
    # If a single value of maximum size is given, set it as the maximum size of all the selections.
  }else if(length(max_size)==1){  
    max_size <- rep(max_size, K)
    # Re-setting maximum sizes that are larger than the sizes of their corresponding predictor pools to avoid breakdown
    if(!all(max_size < pool_size)){
      warning('Some maximum sizes of selections are larger than the sizes of their predictor pools. Re-setting them to the size of the predictor pools.')
      max_size[max_size > pool_size] <- pool_size[max_size > pool_size]
    }
    
  }else{
    # The number of maximum sizes provided must equal to the size of the ensemble
    if(length(max_size)!=K) stop('The length of max_size is not the same as K')
    # Re-setting maximum sizes that are larger than the sizes of their corresponding predictor pools to avoid breakdown
    if(!all(max_size < pool_size)){
      warning('Some maximum size of selections are larger than the size of the predictor pools. Re-setting them to the size of the predictor pools.')
      max_size[max_size > pool_size] <- pool_size[max_size > pool_size]
    }
    
  }
  
  
  # Search Start
  ############### Define a single st2e process ############### 
  one_st2e_itr <- function(current_sel, i){
    # Propose random starts
    if(do_rand_start){
      # If we choose to start from a random predictor set
      # We first uniformly sample a single proposal size, then we sample
      # multiple proposals of the sampled size and choose the best proposal as 
      # the starting set
      best_proposal <- st2_fwd(Y, X, env, new_obj_fun, sel_fun, current_sel, pred_pool[i,], max_size[i],
                               alpha = 1, beta = 1, lambda = 1, a = 1, b = 1, t = 0, fixed_size = T)        
      cur_score <- best_proposal$score
      current_sel <- best_proposal$step
      rand_start <- current_sel
      
    }else{
      cur_score <- new_obj_fun(Y, X, env, current_sel)
      rand_start <- NULL
      
    }
    
    # Choosing maximum allowed consecutive failure of updates till the break of the algorithm,
    # depending on whether we implement switch updates
    max_counter <- ifelse(do_switch, 3, 2)
    counter <- 0
    
    
    while(1){
      # Switch step
      # Calculate the number of predictors that are in and out the model
      n_in <- sum(current_sel==1) 
      n_out <- sum(pred_pool[i,]==1) - n_in
      
      # If there exist predictors to be switched, and the current model size is not larger
      # than the maximum allowed size, do switching as indicated
      if(n_in < max_size[i] & n_in > 0 & n_out > 0 & do_switch){
        best_proposal <- st2_swt(Y, X, env, new_obj_fun, sel_fun, current_sel, pred_pool[i,], max_size[i],       
                                 alpha = alpha, beta = beta, lambda = lambda, a = a, b = b, t = t)
        # Update if the best proposal obtained a lower objective score than the global minimum 
        if(best_proposal$score < cur_score){
          counter <- 0 # Reset the counter for a success update
          cur_score <- best_proposal$score
          current_sel <- best_proposal$step
          
        }else if(best_proposal$score > cur_score){
          counter <- counter + 1
          
        }
        
      }
      # Break the algorithm if reach the maximum allowed consecutive failure of updates
      if(counter >= max_counter) break
      
      # Forward step
      # Calculate the number of predictors that are in the model
      n_in <- sum(current_sel==1) 
      
      # If the current model size is not larger than the maximum allowed size, do forward search
      if(n_in < max_size[i]){
        best_proposal <- st2_fwd(Y, X, env, new_obj_fun, sel_fun, current_sel, pred_pool[i,], max_size[i],       
                                 alpha = alpha, beta = beta, lambda = lambda, a = a, b = b, t = t)
        # Update if the best proposal obtained a lower objective score than the global minimum
        if(best_proposal$score <= cur_score){
          counter <- 0 # Reset the counter for a success update
          cur_score <- best_proposal$score
          current_sel <- best_proposal$step
          
        }else if(best_proposal$score > cur_score){
          counter <- counter + 1
          
        }
        
      }
      # Break the algorithm if reach the maximum allowed consecutive failure of updates
      if(counter >= max_counter) break
      
      # Backward step
      # Calculate the number of predictors that are in the model
      n_in <- sum(current_sel==1) 
      
      # If there exist predictors to be removed, do backward search
      if(n_in > 0){
        best_proposal <- st2_bwd(Y, X, env, new_obj_fun, sel_fun, current_sel, pred_pool[i,], max_size[i],    
                                 alpha = alpha, beta = beta, lambda = lambda, a = a, b = b, t = t)  
        # Update if the best proposal obtained a lower objective score than the global minimum
        if(best_proposal$score < cur_score){
          # If the score of the best blanket is lower then the global minimum, update, otherwise add a counter
          counter <- 0
          cur_score <- best_proposal$score
          current_sel <- best_proposal$step
          
        }else if(best_proposal$score >= cur_score){
          counter <- counter + 1
          
        }
        
      }
      # Break the algorithm if reach the maximum allowed consecutive failure of updates
      if(counter >= max_counter) break
      
    } # This is the end of one ST2 loop
    list(sel = current_sel, score = cur_score,
         model = if(ret_mod) reg_fun(Y, X, env, current_sel[-1]) else NULL,
         start = rand_start)
  }
  
  
  ############### Run st2e without parallel computing ############### 
  if(is.null(ncore)){
    # Retrieved from @https://www.dummies.com/article/technology/programming-web-design/r/how-to-generate-your-own-error-messages-in-r-175112/
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = K,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
    
    for(i in 1:K){
      pb$tick() # Timing progress bar
      res <- one_st2e_itr(ensemble[i,], i)
      ensemble[i,] <- res$sel
      scores[i] <- res$score
      models[[i]] <- res$model
      if(do_rand_start) start_set[i,] <- res$start
    }
    
    
  }else{
    par_method <- match.arg(par_method)
    `%dopar%` <- getFromNamespace("%dopar%", "foreach")
    pkg_export <- if(par_method != 'MC') c(fun_export, 'Matrix') else NULL
    fun_export <- if(par_method != 'MC') c(fun_export, 'rbeta_binom','dbeta_binom','st2_fwd','st2_bwd','st2_swt') else NULL
    
    # Set up cluster
    ############### Run st2e with doParallel, which should work on all platforms ###############
    # SNOW is preferred on a Windows Machine
    if(par_method == 'SNOW'){
      cl<- parallel::makeCluster(ncore)
      doParallel::registerDoParallel(cl)
      
      
    }else if(par_method == 'MC'){
      ############### Run st2e with MC ###############
      # On a non-Windows Machine (e.g. Mac OS, Linux), MC essentially performs the multicore
      # functionality in the Parallel package
      doParallel::registerDoParallel(cores = ncore)
      
      
    }else{
      ############### Run st2e with MPI ############### 
      # MPI should be used when running on a HPC cluster, and
      # should be ran together with command line mpirun in shell scripts.
      if(is.null(doMPI::getDoMpiCluster())){
        cl <- doMPI::startMPIcluster(ncore)
        doMPI::registerDoMPI(cl)
        
      }
      
    }
    
    # Foe each loop
    if(chunk_size == 1){
      res <- foreach::foreach(i=1:K, .verbose = verbose, .packages = pkg_export, .export = fun_export) %dopar% {
        one_st2e_itr(ensemble[i,], i)
        
      }
      
    }else{
      idx_chunk <- split(1:K, ceiling(seq_along(1:K)/chunk_size))
      
      res <- foreach::foreach(i=1:length(idx_chunk), .verbose = verbose, .packages = pkg_export, .export = fun_export) %dopar% {
        out <- list(sel = matrix(nrow = length(idx_chunk[[i]]), ncol = ncol(ensemble)),
                    start = if(do_rand_start) matrix(nrow = length(idx_chunk[[i]]), ncol = ncol(ensemble)) else NULL,
                    score = c(), model = list())
        for(j in 1:length(idx_chunk[[i]])){
          res <- one_st2e_itr(ensemble[idx_chunk[[i]][j],], idx_chunk[[i]][j])
          out$sel[j,] <- res$sel
          out$score[j] <- res$score
          out$model[[j]] <- res$model
          if(do_rand_start) out$start[j,] <- res$start
          
        }
        out
        
      }
      
    }
    
    ensemble <- do.call(rbind, lapply(res, function(x) x$sel))
    colnames(ensemble) <- pred_names
    scores <- do.call(c, lapply(res, function(x) x$score))
    models <-  if(ret_mod) do.call(c, lapply(res, function(x) x$model)) else NULL
    if(do_rand_start) start_set <- do.call(rbind, lapply(res, function(x) x$start)) %>% Matrix::Matrix()
    if(par_method == 'SNOW') parallel::stopCluster(cl)
    
  }
  
  ############### Prepare the output ############### 
  ensemble <- Matrix::Matrix(ensemble)
  out <- list(ensemble = ensemble, obj_scores = scores, models = models, 
              psd_pred = psd_pred, 
              pred_pool = pred_pool, start_set = start_set,
              imp_scores = NULL)
  class(out) <- 'st2e'
  
  
  ############### Calculate importance scores of predictors ############### 
  if(ret_imp){
    # A list containing all the default controls for calculating importance scores
    default <- list(prune = T, prefilt_scores = NULL, scale_psd_imp = F, per_pred_evl = F, B = 5000)
    
    # Check whether all the controls provided for calculating importance scores are valid
    ctrl_name <- names(calc_imp_ctrl)
    default_name <- names(default)
    unknown_name <- ctrl_name[!ctrl_name %in% default_name]
    
    # Warning if there is an unrecognized control
    if(length(unknown_name) >= 1){
      warning("Unrecognized names in the cal_imp control: ", paste(unknown_name, collapse=", "), 
              ", These unrecognized controls are discarded.")
      calc_imp_ctrl <- calc_imp_ctrl[!ctrl_name %in% unknown_name]
      ctrl_name <- names(calc_imp_ctrl)
      
    }
    
    # Update controls
    default[ctrl_name] <- calc_imp_ctrl
    
    # Calculate importance
    out <- do.call('calc_imp', c(st2e_obj = list(out), default))
    
  }
  
  out
  
}


# The function for making forward selections in ST2
st2_fwd <- function(Y, X, env, obj_fun, sel_fun, pred_in, pred_pool,
                    max_size, alpha, beta, lambda, a, b, t, fixed_size = F){
  
  # The size of the predictor pool
  pool_size <- sum(pred_pool)
  
  # Obtain the number and the indexes of predictors that are not included
  idx_out <- which(pred_pool - pred_in ==1)
  n_out <- length(idx_out)
  n_in <-  sum(pred_in)
  
  # Calculate the maximum number of predictors allowed to be added in this step
  # Refer to Supplementary Method 7.1.4 eq. (6) for details
  n_add_alw <- max_size - n_in
  
  # Calculate the maximum number of predictors to be added in this step
  nu <- ceiling(lambda*(n_add_alw- 1))
  if(fixed_size){
    # If sample proposals of a same size, we first sample a size 
    delta <- rbeta_binom(1, nu, alpha, beta) + 1
    # Calculate the number of proposals to sample
    # Refer to Supplementary Method 7.1.4 eq. (8) for details
    eta <- ceiling(a * n_out ^ t +
                     b * lchoose(n_out, delta))
    # Replicate delta to make it an vector of eta
    delta <- rep(delta, eta)
    
  }else{
    # Calculate the number of proposals to sample
    # Refer to Supplementary Method 7.1.4 eq. (8) for details
    eta <- ceiling(a * n_out ^ t +
                     b * sum(
                       dbeta_binom(0:nu, nu, shape1 = alpha, shape2 = beta) * lchoose(n_out, (0:nu + 1))
                     )
    )
    delta <- rbeta_binom(eta, nu, alpha, beta) + 1 # Sample eta numbers of proposal sizes
    
  }
  
  # Sample eta numbers of proposals
  proposals <- matrix(rep(pred_in, eta), nrow = eta, byrow = T) 
  if(length(idx_out)==1){
    # If there is only one predictor outside of the model, add it to the model
    proposals[,idx_out] <- 1
    
  }else{
    # At each iteration, sample a proposal of size delta[i] from predictors outside of the model
    for(i in 1:eta) proposals[i,sample(idx_out,size = delta[i])] <- 1
    
  }
  
  # Remove duplicated proposals
  proposals <- proposals[!duplicated(proposals),,drop = FALSE] 
  
  eta <- nrow(proposals)
  scores <- c()
  
  # Calculate objective scores for each proposals
  for(i in 1:eta){
    scores[i] <- obj_fun(Y, X, env, proposals[i,]) 
    
  } 
  
  # Select the best proposal
  best <- sel_fun(scores) 
  list(score = scores[best], step = proposals[best,]) 
  
}


# The function for making backward selections in ST2
st2_bwd <- function(Y, X, env, obj_fun, sel_fun, pred_in, pred_pool,
                    max_size, alpha, beta, lambda, a, b, t){
  
  # The size of the predictor pool
  pool_size <- sum(pred_pool)
  
  # Obtain the number and the indexes of predictors that are not included
  n_in <- sum(pred_in)
  idx_in <- which(pred_in == 1)
  
  # Calculate the number of predictors that must be removed in this step
  # Refer to Supplementary Method 7.1.4 eq. (7) for details
  n_must_rmv <- max(n_in - max_size, 1)
  
  # Calculate the maximum number of predictors to be removed in this step
  nu <- ceiling(lambda*(n_in - n_must_rmv))
  
  # Calculate the number of proposals to sample
  # Refer to Supplementary Method 7.1.4 eq. (9) for details
  eta <-  ceiling(a * n_in ^ t + 
                    b * sum(
                      dbeta_binom(0:nu, nu, shape1 = alpha, shape2 = beta) * lchoose(n_in, (0:nu + n_must_rmv))
                    )
  )
  
  delta <- rbeta_binom(eta, nu, alpha, beta) + n_must_rmv # Sample eta numbers of proposal sizes
  
  # Sample eta numbers of proposals
  proposals <- matrix(rep(pred_in, eta), nrow = eta, byrow = T)
  if(length(idx_in)==1){
    # If there is only one predictor in the model, remove it from the model
    proposals[,idx_in] <- 0
    
  }else{
    # At each iteration, sample a proposal of size delta[i] from predictors in the model
    for(i in 1:eta) proposals[i,sample(idx_in, delta[i])] <- 0
    
  }
  
  # Remove duplicated proposals
  proposals <- proposals[!duplicated(proposals),,drop = FALSE] 
  
  eta <- nrow(proposals)
  scores <- c()
  
  # Calculate the objective score for each proposal
  for(i in 1:eta){
    scores[i] <- obj_fun(Y, X, env, proposals[i,])
    
  } 
  
  # Select the best proposal
  best <- sel_fun(scores) 
  list(score = scores[best], step = proposals[best,]) 
  
}


# The function for making switching of predictors in ST2
st2_swt <- function(Y, X, env, obj_fun, sel_fun, pred_in, pred_pool,
                    max_size, alpha, beta, lambda, a, b, t){
  
  # The size of the predictor pool
  pool_size <- sum(pred_pool)
  
  # Obtain the number and the indexes of predictors that are not included
  idx_out <- which(pred_pool - pred_in == 1)
  n_out <- length(idx_out)
  
  # Obtain the number and the indexes of predictors that are included
  idx_in <- which(pred_in == 1)
  n_in <-  length(idx_in)
  
  # Calculate the maximum number of predictors allowed to be switched in this step
  n_swt_alw <- min(n_out, n_in)
  
  # Calculate the maximum number of predictors to be switched in this step
  nu <- ceiling(lambda*(n_swt_alw - 1))
  
  # Calculate the number of proposals to sample
  eta <- ceiling(a * n_in ^ t +
                   b * sum(
                     dbeta_binom(0:nu, nu, shape1 = alpha, shape2 = beta) * 
                       (lchoose(n_in, (0:nu + 1)) + lchoose(n_out, (0:nu + 1)))
                   )
  )
  
  delta <- rbeta_binom(eta, nu, alpha, beta) + 1 # Sample eta numbers of proposal sizes
  
  # Sample eta numbers of proposals to remove from the model
  proposals <- matrix(rep(pred_in, eta),nrow = eta, byrow = T) 
  if(length(idx_in)==1){
    # If there is only one predictor in the model, remove it from the model
    proposals[,idx_in] <- 0
    
  }else{
    # At each iteration, sample a proposal of size delta[i] from predictors in the model
    for(i in 1:eta) proposals[i,sample(idx_in,size = delta[i])] <- 0
    
  }
  # Sample eta numbers of proposals to add into the model
  if(length(idx_out)==1){
    # If there is only one predictor outside of the model, add it to the model
    proposals[,idx_out] <- 1
    
  }else{
    # At each iteration, sample a proposal of size delta[i] from predictors outside the model
    for(i in 1:eta) proposals[i,sample(idx_out,size = delta[i])] <- 1
    
  }
  
  # Remove duplicated proposals
  proposals <- proposals[!duplicated(proposals),,drop = FALSE]
  
  eta <- nrow(proposals)
  scores <- c()
  
  # Calculate the objective score for each proposal
  for(i in 1:eta){
    scores[i] <- obj_fun(Y, X, env, proposals[i,]) 
    
  } 
  
  # Select the best proposal
  best <- sel_fun(scores)
  list(score = scores[best], step = proposals[best,]) 
  
}


# StableMate
stablemate <- function(Y, X, env, K = 100, 
                       alpha = 1, beta = 5, lambda = 0.5, a = 1, b = 1, t = 0, max_size = NULL, 
                       pred_st2_ctrl = NULL, stab_st2_ctrl = NULL,
                       pred_calc_imp_ctrl = NULL, stab_calc_imp_ctrl = NULL,
                       do_switch = FALSE, ret_mod = TRUE, ret_imp = TRUE, ncore = NULL,
                       par_method = c('SNOW','MC','MPI'), verbose = TRUE, 
                       chunk_size = ceiling(K/ncore), fun_export = NULL, pkg_export = NULL,
                       drop_pred = TRUE){
  ############### Check whether the ST2 stability selection is based on the ST2 predictivity selection############### 
  # If the selections of stable predictors are based on the prediction ensemble built
  # by the first run of ST2E, then calculate importance scores of predictors according to 
  # Supplementary Methods 7.1.3 eq.(7,8,9) where the importance scores measuring predictivity and
  # stability are matched in the sense that the importance scores measuring stability are 
  # calculated derived from the importance scores measuring predictivity
  match_sel <- TRUE
  if(!is.null(stab_calc_imp_ctrl$pred_pool)){
    if(!is.vector(stab_calc_imp_ctrl$pred_pool)){
      match_sel <- FALSE
    }else if(length(stab_calc_imp_ctrl$pred_pool)!=1){
      match_sel <- FALSE
    }else if(stab_calc_imp_ctrl$pred_pool != 'prediction_ensemble'){
      match_sel <- FALSE
    }
  }
  
  if(!is.null(stab_calc_imp_ctrl$start_set)){
    if(!is.vector(stab_calc_imp_ctrl$start_set)){
      match_sel <- FALSE
    }else if(length(stab_calc_imp_ctrl$start_set)!=1){
      match_sel <- FALSE
    }else if(!stab_calc_imp_ctrl$start_set %in% c('prediction_ensemble','random_start')){
      match_sel <- FALSE
    }
  }
  
  if(!match_sel & drop_pred){
    warning('The stability selection is not based on the predictivity selection. Reset drop_pred to FALSE.')
    drop_pred <- FALSE
  }
  
  
  ############### Initialize the ST2 controls for selecting predictive sets ############### 
  # A list containing all the default ST2 controls for predictivity selections
  defaults <- list(obj_fun = obj_bic, sel_fun = sel_min, 
                   reg_fun = reg_ols, pred_pool = NULL, start_set = 'random_start',
                   alpha = alpha, beta = beta, lambda = lambda,
                   a = a, b = b, t = t, max_size = max_size, 
                   do_switch = do_switch, ret_mod = ret_mod, ret_imp = FALSE, K = K, calc_imp_ctrl = NULL,
                   ncore = ncore, par_method = par_method, verbose = verbose,
                   chunk_size = chunk_size, fun_export = fun_export, pkg_export = pkg_export,
                   drop_pred = drop_pred)
  
  # Check whether all the controls provided for selecting predictive sets selection are valid
  ctrl_names <- names(pred_st2_ctrl)
  default_names <- names(defaults)
  unknown_names <- ctrl_names[!ctrl_names %in% default_names]
  
  # Set global controls. These controls cannot be adjusted by pred_st2_ctrl
  global_ctrls <- c('K', 'ret_imp', 'calc_imp_ctrl','drop_pred')
  if(any(ctrl_names %in% global_ctrls)){
    warning("Discard ", paste(ctrl_names[ctrl_names %in% global_ctrls], collapse=", "),
            " in pred_st2_ctrl, use global control instead.")
    pred_st2_ctrl <- pred_st2_ctrl[!ctrl_names %in% global_ctrls]
    ctrl_names <- names(pred_st2_ctrl)
    
  }
  
  # Warning if there is an unrecognized control
  if(length(unknown_names) >= 1){
    warning("Unrecognized names in the st2e control (predictivity): ", paste(unknown_names, collapse=", "), 
            ", These unrecognized controls are discarded.")
    pred_st2_ctrl <- pred_st2_ctrl[!ctrl_names %in% unknown_names]
    ctrl_names <- names(pred_st2_ctrl)
    
  } 	
  
  # Update the controls
  defaults[ctrl_names] <- pred_st2_ctrl
  defaults$chunk_size <- ceiling(defaults$K/defaults$ncore)
  
  ############### Build prediction ensemble ############### 
  mod_pred <- do.call('st2e', 
                      c(list(Y = Y, X = X, env = env), 
                        defaults))
  
  
  ############### Initialize the ST2 controls for selecting stable sets ############### 
  # A list containing all the default ST2 controls for stability selections
  defaults <- list(obj_fun = obj_pss, sel_fun = sel_min,
                   reg_fun = reg_ols, pred_pool = 'prediction_ensemble', start_set = 'random_start',
                   alpha = alpha, beta = beta, lambda = lambda,
                   a = a, b = b, t = t, max_size = max_size,
                   do_switch = do_switch, ret_mod = ret_mod, ret_imp = FALSE, K = K, calc_imp_ctrl = NULL,
                   ncore = ncore, par_method = par_method, verbose = verbose,
                   chunk_size = chunk_size, fun_export = fun_export, pkg_export = pkg_export,
                   drop_pred = FALSE)
  # Check whether all the controls provided for selecting stable sets are valid
  ctrl_names <- names(stab_st2_ctrl)
  default_names <- names(defaults)
  unknown_names <- ctrl_names[!ctrl_names %in% default_names]
  
  # Set global controls. These controls cannot be adjusted by stab_st2_ctrl
  global_ctrls <- c('K', 'ret_imp', 'calc_imp_ctrl','drop_pred')
  if(any(ctrl_names %in% global_ctrls)){
    warning("Discard ", paste(ctrl_names[ctrl_names %in% global_ctrls], collapse=", "),
            " in stab_st2_ctrl, use global control instead.")
    stab_st2_ctrl <- stab_st2_ctrl[!ctrl_names %in% global_ctrls]
    ctrl_names <- names(stab_st2_ctrl)
    
  }
  
  # Warning if there is an unrecognized control
  if(length(unknown_names) >= 1){
    warning("Unrecognized names in the st2e control (stability): ", paste(unknown_names, collapse=", "), 
            ", These unrecognized controls are discarded.")
    stab_st2_ctrl <- stab_st2_ctrl[!ctrl_names %in% unknown_names]
    ctrl_names <- names(stab_st2_ctrl)
    
  }
  
  # Update the controls
  defaults[ctrl_names] <- stab_st2_ctrl
  defaults$chunk_size <- ceiling(defaults$K/defaults$ncore)
  
  ############### Build stable ensemble ############### 
  message('Now search for stable sets')
  
  if(match_sel){
    # Initialize the predictor pools for predictivity selection
    defaults$pred_pool <- mod_pred$ensemble
    defaults$pred_pool[,1] <- 1
    
    # Initialize the starting sets for stability selection
    if(defaults$start_set == 'prediction_ensemble'){
      defaults$start_set <- mod_pred$ensemble
      defaults$start_set[,1] <- 1
    }
  }
  
  if(drop_pred){
    pred_in_pool <- colnames(mod_pred$pred_pool)
    X <- X[,pred_in_pool[-1]]
  }
  
  # Run ST2E
  mod_stab <- do.call('st2e', 
                      c(list(Y = Y, X = X, env = env), 
                        defaults))
  
  
  if(ret_imp){
    
    if(match_sel){
      ############### Calculate importance scores ############### 
      message('The stable ensemble is built based on the prediction ensemble: the importance scores measuring predictivity and stability should be calculated together and matched. Use the pred_calc_imp_ctrl as the control for cal_imp function.')
      # If predictor pools provided are the results of StableMate pre-filtering based on Lasso, scale the
      # importance score of the pseudo-predictor
      if(!is.null(attributes(pred_st2_ctrl$pred_pool)$stbm_lasso) & is.null(pred_calc_imp_ctrl$scale_psd_imp)){
        pred_calc_imp_ctrl$scale_psd_imp <- TRUE
        
      }
      
      # A list containing all the default controls for calculating importance scores
      defaults <- list(prune = T, prefilt_scores = NULL, scale_psd_imp = F, per_pred_evl = F, B = 5000)
      
      # Check whether all the controls provided for calculating importance scores are valid
      ctrl_names <- names(pred_calc_imp_ctrl)
      default_names <- names(defaults)
      unknown_names <- ctrl_names[!ctrl_names %in% default_names]
      
      # Warning if there is an unrecognized control
      if(length(unknown_names) >= 1){
        warning("Unrecognized names in the cal_imp control: ", paste(unknown_names, collapse=", "), 
                ", These unrecognized controls are discarded.")
        pred_calc_imp_ctrl <- pred_calc_imp_ctrl[!ctrl_names %in% unknown_names]
        ctrl_names <- names(pred_calc_imp_ctrl)
        
      } 	
      
      # Update the controls
      defaults[ctrl_names] <- pred_calc_imp_ctrl
      
      # Calculate importance
      out <- do.call('calc_imp', c(st2e_obj = list(mod_stab), st2e_obj2 = list(mod_pred), defaults))
      names(out) <- c('stable_ensemble','prediction_ensemble')
      out$matched <- match_sel
      
      # Calculate the importance scores measuring predictivity and stability separately
    }else{
      warning("The stable ensemble is not built based on the prediction ensemble. Calculate the importance scores measuring predictivity and stability seperately.")
      
      ############### Calculate the importance scores measuring predictivity of predictors ############### 
      # If predictor pools provided are the results of StableMate pre-filtering based on Lasso, scale the
      # importance score of the pseudo-predictor
      if(!is.null(attributes(pred_st2_ctrl$pred_pool)$stbm_lasso) & is.null(pred_calc_imp_ctrl$scale_psd_imp)){
        pred_calc_imp_ctrl$scale_psd_imp <- TRUE
        
      }
      
      # A list containing all the default controls for calculating importance scores measuring predictivity
      defaults <- list(prune = T, prefilt_scores = NULL, scale_psd_imp = F, per_pred_evl = F, B = 5000)
      
      # Check whether all the controls provided for calculating importance scores are valid
      ctrl_names <- names(pred_calc_imp_ctrl)
      default_names <- names(default_names)
      unknown_names <- ctrl_names[!ctrl_names %in% default_names]
      
      # Warning if there is an unrecognized control
      if(length(unknown_names) >= 1){
        warning("Unrecognized names in the cal_imp control (predictivity): ", paste(unknown_names, collapse=", "), 
                ", These unrecognized controls are discarded.")
        pred_calc_imp_ctrl <- pred_calc_imp_ctrl[!ctrl_names %in% unknown_names]
        ctrl_names <- names(pred_calc_imp_ctrl)
        
      }
      
      # Update controls
      defaults[ctrl_names] <- pred_calc_imp_ctrl
      
      # Calculate importance
      mod_pred <- do.call('calc_imp', c(st2e_obj = list(mod_pred), defaults))
      
      
      ############### Calculate the importance scores of measuring stability of predictors ############### 
      # If predictor pools provided are the results of StableMate pre-filtering based on Lasso, scale the
      # importance score of the pseudo-predictor
      if(!is.null(attributes(stab_st2_ctrl$pred_pool)$stbm_lasso) & is.null(stab_calc_imp_ctrl$scale_psd_imp)){
        stab_calc_imp_ctrl$scale_psd_imp <- TRUE
        
      }
      
      # A list containing all the default controls for calculating the importance scores measuring stability of predictors
      defaults <- list(prune = T, prefilt_scores = NULL, scale_psd_imp = F, per_pred_evl = F, B = 5000)
      
      # Check whether all the controls provided for calculating importance scores are valid
      ctrl_names <- names(stab_calc_imp_ctrl)
      default_names <- names(defaults)
      unknown_names <- ctrl_names[!ctrl_names %in% default_names]
      
      # Warning if there is an unrecognized control
      if(length(unknown_names) >= 1){
        warning("Unrecognized names in cal_imp control (stability): ", paste(unknown_names, collapse=", "), 
                ", These unrecognized controls are discarded.")
        stab_calc_imp_ctrl <- stab_calc_imp_ctrl[!ctrl_names %in% unknown_names]
        ctrl_names <- names(stab_calc_imp_ctrl)
        
      } 	
      
      # Update the controls
      defaults[ctrl_names] <- stab_calc_imp_ctrl
      
      # Calculate importance
      mod_stab <- do.call('calc_imp', c(st2e_obj = list(mod_stab), defaults))
      
      out <- list(stable_ensemble = mod_stab, prediction_ensemble = mod_pred, matched = match_sel)
      
    }
    
  }else{
    out <- list(stable_ensemble = mod_stab, prediction_ensemble = mod_pred, matched = match_sel)
    
  }
  
  # Return
  class(out) <- 'stablemate'
  out
  
}


# Calculate importance scores
calc_imp <- function(st2e_obj, st2e_obj2 = NULL, prune = T, prefilt_scores = NULL, scale_psd_imp = F,
                     per_pred_evl = F, B = 1000){
  
  if(class(st2e_obj) != 'st2e') stop("The st2e_obj provided is not a st2e object")
  
  K <- nrow(st2e_obj$ensemble) # Size of the ensemble
  P <- ncol(st2e_obj$ensemble) # The number of predictors
  pred_names <- colnames(st2e_obj$ensemble) # Predictor names
  
  # If st2e_obj is obtained by running ST2E again based on the st2e_obj2 selection results (and if st2e_obj2 is given),
  # check whether the predictor pool of st2e_obj matches with the ensemble of st2e_obj2
  matched_psd <- FALSE
  if(!is.null(st2e_obj2)){
    if(class(st2e_obj2) != 'st2e'){
      stop("The st2e_obj2 provided is not a st2e object, ignore st2e_obj2.")
      
    }else if(!all(st2e_obj$pred_pool[,-1] == st2e_obj2$ensemble[,-1])){
      warning("The selections of st2e_obj are not based on the ensemble of st2e_obj2, ignore st2e_obj2.")
      st2e_obj2 <- NULL
      
    }else if(!all(st2e_obj$pred_pool[,1]==1) & !all(st2e_obj$pred_pool[,1] == st2e_obj2$ensemble[,1])){
      warning("The pseudo-predictor should be either included in all the predictor pools of st2e_obj, or match with the selections of st2e_obj2, ignore st2e_obj2.")
      st2e_obj2 <- NULL
      
    }else if(all(st2e_obj$pred_pool[,1] == st2e_obj2$ensemble[,1])){
      matched_psd <- TRUE
      
    }
    
  } 
  
  
  ################ Calculate weights of selections  ################ 
  # Calculate the weights of st2e_obj selections if prune the ensemble (In stablemate, its the ranking quantiles of pss scores)
  if(prune) q_sel <- rank(-st2e_obj$obj_scores, ties.method = 'average')/K else q_sel <- rep(1,K)
  # Calculate the weights of st2e_obj2 selections if prune the ensemble (In stablemate, its the ranking quantiles of bic scores)
  if(prune & !is.null(st2e_obj2)) q_sel2 <- rank(-st2e_obj2$obj_scores, ties.method = 'average')/K else q_sel2 <- rep(1,K)
  # Calculate the weights of prior selections if prune the ensemble (In stablemate, its the ranking quantiles of the objective scores of the pre-filtering method applied)
  if(prune & !is.null(prefilt_scores)) q_prefilt <- rank(-prefilt_scores, ties.method = 'average')/K else q_prefilt <-  rep(1,K)
  
  # Weight selections by whether they are both important for the pre-filtering and st2e_obj2 (inspired by probability t-norm)
  q_sel2 <- q_sel2*q_prefilt
  # Weight selections by whether they are both important for the pre-filtering, st2e_obj1 and st2e_obj2 (inspired by probability t-norm)
  q_sel <- q_sel*q_sel2
  
  # Constants for normalizing importance scores so that they range between [0,1]
  q_sel_sum <- sum(q_sel)
  q_sel_sum2 <- sum(q_sel2)
  q_prefilt_sum <- sum(q_prefilt)
  
  
  ################ Calculate the importance scores of the predictors selected by st2e_obj  ################ 
  # Calculate the joint importance scores of the predictors selected by st2e_obj
  # (In stablemate, its a measure of stability and predictivity of predictors,
  # refer to Supplementary Methods 7.1.3 eq.(7))
  sel_imp <- colSums(st2e_obj$ensemble*q_sel)/q_sel_sum
  # Set up the matrix for storing the bootstrapped joint importance scores in st2e_obj
  bt_sel_imp <- matrix(0, nrow=B, ncol= P)
  colnames(bt_sel_imp) <- pred_names
  
  # Calculate prior importance scores of the predictors selected by st2e_obj,
  # which are also the joint importance scores of the predictors selected by st2e_obj2
  # (In stablemate, its a measure of predictivity of predictors,
  # refer to Supplementary Methods 7.1.3 eq.(8))
  prior_imp <- colSums(st2e_obj$pred_pool*q_sel2)/q_sel_sum2
  
  # Calculate the conditional importance scores of the predictors selected by st2e_obj
  # (In stablemate, its a measure of stability of predictors,
  # refer to Supplementary Methods 7.1.3 eq.(9))
  condi_imp <- pmin(1, sel_imp/prior_imp)
  names(condi_imp) <- pred_names
  # Set up the matrix for storing the bootstrapped conditional importance scores in st2e_obj
  bt_condi_imp <- matrix(0,nrow=B, ncol=P)
  colnames(bt_condi_imp) <- pred_names
  
  # Set up the predictor-wise evaluation of significance of conditional importance scores,
  # in which we generate predictor-specific benchmark for each predictor based on the 
  # pseudo-predictor
  
  if(per_pred_evl){
    if(all(st2e_obj$pred_pool[,1] == 1)){
      # Build a matrix that replicates P selections results of the pseudo-predictor over K repeated selections
      mask_mat <- matrix(rep(st2e_obj$ensemble[,1],P),ncol = P)
      # Multiply it by the predictor pool, we get for each true predictor, a selection result of the pseudo-predictor within
      # the pool of the true predictor. We are going to compare the selection of the pseudo-predictor 
      # and the true predictor within the pool of the true predictor
      mask_mat <- mask_mat * st2e_obj$pred_pool
      colnames(mask_mat) <- pred_names
      # Create a matrix for storing the bootstrapped conditional importance scores of the pseudo-predictor
      benchmark <- matrix(0,nrow=B, ncol=P)
      colnames(benchmark) <- pred_names
      
    }else{
      warning("The pseudo-predictor is not included in all the predictor pools of st2e_obj, set per_pred_evl to FALSE.")
      per_pred_evl <- FALSE
    }
    
  }
  
  ################ Calculate the importance scores of the predictors selected by st2e_obj2  ################ 
  if(!is.null(st2e_obj2)){
    # Calculate the joint importance scores of the predictors selected by st2e_obj2
    sel_imp2 <- prior_imp
    
    if(!matched_psd){
      # Adjust the joint importance score of the pseudo-predictor in st2e_obj2
      sel_imp2[1] <- sum(st2e_obj2$ensemble[,1]*q_sel2)/q_sel_sum2
      # Adjust the joint importance score of the pseudo-predictor in st2e_obj
      sel_imp[1] <- sum(st2e_obj2$ensemble[,1]*st2e_obj$ensemble[,1]*q_sel)/q_sel_sum
      
    }
    
    # Calculate importance scores of predictors in pre-filtering based
    # on the predictor pools of st2e_obj2 if st2e_obj2 is given
    prefilt_imp <- colSums(st2e_obj2$pred_pool*q_prefilt)/q_prefilt_sum 
    
    # Set up the matrix for storing the bootstrapped joint importance scores in st2e_obj2
    bt_sel_imp2 <- matrix(0, nrow=B, ncol= P)
    colnames(bt_sel_imp2) <- pred_names
    
    # Calculate the conditional importance scores of the predictors selected by st2e_obj2
    condi_imp2 <- pmin(1, sel_imp2/prefilt_imp)
    names(condi_imp2) <- pred_names
    # Set up the matrix for storing the bootstrapped conditional importance scores in st2e_obj2
    bt_condi_imp2 <- matrix(0, nrow=B, ncol= P)
    colnames(bt_condi_imp2) <- pred_names
    
    # Set up predictor-wise evaluation of the significance of the conditional importance scores the predictors selected by st2e_obj2,
    # similar as what has been done for st2e_obj
    if(per_pred_evl){
      if(all(st2e_obj2$pred_pool[,1] == 1)){
        mask_mat2 <- matrix(rep(st2e_obj2$ensemble[,1],P),ncol = P)
        mask_mat2 <- mask_mat2 * st2e_obj2$pred_pool
        colnames(mask_mat2) <- pred_names
        benchmark2 <- matrix(0,nrow=B, ncol=P)
        colnames(benchmark2) <- pred_names
        
      }else{
        warning("The pseudo-predictor is not included in all the predictor pools of st2e_obj2, set per_pred_evl to FALSE.")
        per_pred_evl <- FALSE
      }
      
    }
    
  }
  
  
  ################ Scale the importance scores of the pseudo-predictor ################ 
  # Scale the prior and the joint importance scores of the pseudo-predictor if the pre-filtering process 
  # cannot apply proper selections on the pseudo-predictor. Assume that all the selections in the ensemble
  # were applied with the same stochastic pre-filtering procedure, and assume that the pseudo-predictor is
  # included in all the predictor pools
  if(scale_psd_imp){
    if(!is.null(st2e_obj2)){
      # First calculate sizes of predictor pools
      pool_sizes <- rowSums(st2e_obj2$pred_pool[,-1]) 
      # Calculate how many predictors were retained in average across the ensemble after pre-filtering
      avg_pool_size <- floor(mean(pool_sizes))
      # Check whether the pseudo-predictor is included in all the predictor pools
      if(!all(st2e_obj2$pred_pool[,1] == 1)) warning("The pseudo-predictor is not included in all the predictor pools of st2e_obj2. The calculation of joint importance may be invalid.")
      
      # The importance score of the pseudo-predictor will be multiplied by the avg_pool_size(th) largest importance score
      # in pre-filtering 
      s <- sort(prefilt_imp[-1], decreasing = T)[avg_pool_size]
      prefilt_imp[1] <- prefilt_imp[1] * s
      sel_imp2[1] <- sel_imp2[1] * s
      
    }else{
      if(!all(st2e_obj$pred_pool[,1] == 1)) warning("The pseudo-predictor is not included in all the predictor pools of st2e_obj. The calculation of joint importance may be invalid.")
      pool_sizes <-rowSums(st2e_obj$pred_pool[,-1])
      avg_pool_size <- floor(mean(pool_sizes))
      s <- sort(prior_imp[-1], decreasing = T)[avg_pool_size]
      
    }
    
    sel_imp[1] <- sel_imp[1] * s
    
  }
  
  
  ################ Bootstrap selections  ################ 
  for(b in 1:B){
    k <- sample(1:K, K,replace = T) # Sample K indexes
    
    # Calculate weights of bootstrapped selections (Same as above except for bootstrapping)
    if(prune) q_sel <- rank(-st2e_obj$obj_scores[k], ties.method = 'average')/K else q_sel <- rep(1,K)
    if(prune & !is.null(st2e_obj2)) q_sel2 <- rank(-st2e_obj2$obj_scores[k], ties.method = 'average')/K else q_sel2 <- rep(1,K)
    if(prune & !is.null(prefilt_scores)) q_prefilt <- rank(-prefilt_scores[k], ties.method = 'average')/K else q_prefilt <-  rep(1,K)
    
    q_sel2 <- q_sel2*q_prefilt
    q_sel <- q_sel*q_sel2
    q_sel_sum <- sum(q_sel)
    q_sel_sum2 <- sum(q_sel2)
    
    # Bootstrapped joint importance scores
    bt_sel_imp[b,] <- colSums(st2e_obj$ensemble[k,]*q_sel)/q_sel_sum
    # Bootstrapped prior importance scores
    bt_prior_imp <- colSums(st2e_obj$pred_pool[k,]*q_sel2)/q_sel_sum2
    # Bootstrapped conditional importance scores
    bt_condi_imp[b,] <-  bt_sel_imp[b,]/bt_prior_imp
    # Create a matrix storing the bootstrapped conditional importance scores of the pseudo-predictor
    # calculated with respect to the predictor pool of each true predictor
    if(per_pred_evl) benchmark[b,] <- colSums(mask_mat[k,]*q_sel)/q_sel_sum/bt_prior_imp
    
    if(!is.null(st2e_obj2)){
      # Bootstrapped joint importance scores in st2e_obj2
      bt_sel_imp2[b,] <- bt_prior_imp
      if(!matched_psd){
        # Adjust the bootstrapped joint importance score of the pseudo-predictor in st2e_obj2
        bt_sel_imp[b,1] <- sum(st2e_obj2$ensemble[k,1]*st2e_obj$ensemble[k,1]*q_sel)/q_sel_sum
        # Adjust the bootstrapped joint importance score of the pseudo-predictor in st2e_obj2
        bt_sel_imp2[b,1]  <- sum(st2e_obj2$ensemble[k,1]*q_sel2)/q_sel_sum2
        
      }
      
      # Bootstrapped importance scores of the pre-filtered sets
      bt_prefilt_imp <- colSums(st2e_obj2$pred_pool[k,]*q_prefilt)/q_prefilt_sum
      # Bootstrapped conditional importance scores in st2e_obj2
      bt_condi_imp2[b,] <-  bt_sel_imp2[b,]/bt_prefilt_imp
      # Bootstrapped conditional importance scores of the pseudo-predictor
      if(per_pred_evl) benchmark2[b,] <- colSums(mask_mat2[k,]*q_sel2)/q_sel_sum2/bt_prefilt_imp
      
    }
    
    # Scale the importance score of the pseudo-predictor if the pre-filtering procedure
    # cannot filter the pseudo-predictor
    if(scale_psd_imp){
      avg_pool_size <-  floor(mean(pool_sizes[k]))
      if(!is.null(st2e_obj2)){
        s <- sort(bt_prefilt_imp[-1], decreasing = T)[avg_pool_size]
        bt_sel_imp2[b,1] <- bt_sel_imp2[b,1] * s
      }else{
        s <- sort(bt_prior_imp[-1], decreasing = T)[avg_pool_size]
      }
      bt_sel_imp[b,1] <- bt_sel_imp[b,1] * s
      
    } 
    
  }
  
  
  ################ Calculate significance of joint importance  ################ 
  # For a true predictor, calculate how often it obtains higher (lower) joint importance scores 
  # than the pseudo-predictor over the bootstrap iteration
  bt_sel_imp <- Matrix::Matrix(bt_sel_imp)
  significance <- colMeans(bt_sel_imp > bt_sel_imp[,1]) 
  lt_significance <- colMeans(bt_sel_imp < bt_sel_imp[,1])
  
  joint <- list(importance = sel_imp, bt_imp = bt_sel_imp, 
                significance = significance, lt_significance = lt_significance)
  
  
  ################ Calculate significance of conditional importance  ################
  # Clean up NAs in conditional importance scores that are due to zeros in prior
  # importance scores
  condi_imp[is.na(condi_imp)] <- 0
  bt_condi_imp[is.na(bt_condi_imp)] <- 0
  bt_condi_imp <- Matrix::Matrix(bt_condi_imp)
  
  # For a true predictor, calculate how often it obtains higher (lower) conditional importance scores
  # than the pseudo-predictor over the bootstrap iteration
  if(per_pred_evl){
    benchmark[is.na(benchmark)] <- 0
    benchmark <- Matrix::Matrix(benchmark)
    significance <- colMeans(bt_condi_imp > benchmark) 
    lt_significance <- colMeans(bt_condi_imp < benchmark)
    
  }else{
    benchmark <- NULL 
    significance <- colMeans(bt_condi_imp > bt_condi_imp[,1])
    lt_significance <- colMeans(bt_condi_imp < bt_condi_imp[,1])
    
  }
  
  conditional <- list(importance = condi_imp, bt_imp = bt_condi_imp,
                      benchmark = benchmark, significance = significance, lt_significance = lt_significance)
  
  # Store results in the imp_score slot of the ST2E object
  st2e_obj$imp_scores <- list(joint = joint, conditional = conditional)
  
  ################ Similar calculation of significance in st2e_obj2 ################
  if(!is.null(st2e_obj2)){
    # Joint
    bt_sel_imp2 <- Matrix::Matrix(bt_sel_imp2)
    significance <- colMeans(bt_sel_imp2 > bt_sel_imp2[,1]) 
    lt_significance <- colMeans(bt_sel_imp2 < bt_sel_imp2[,1])
    
    joint2 <- list(importance = sel_imp2, bt_imp = bt_sel_imp2, 
                   significance = significance, lt_significance = lt_significance)
    
    # Conditional
    condi_imp2[is.na(condi_imp2)] <- 0
    bt_condi_imp2[is.na(bt_condi_imp2)] <- 0
    bt_condi_imp2 <- Matrix::Matrix(bt_condi_imp2)
    if(per_pred_evl){
      benchmark2[is.na(benchmark2)] <- 0
      benchmark2 <- Matrix::Matrix(benchmark2)
      significance <- colMeans(bt_condi_imp2 > benchmark2) 
      lt_significance <- colMeans(bt_condi_imp2 < benchmark2)
      
    }else{
      benchmark2 <- NULL 
      significance <- colMeans(bt_condi_imp2 > bt_condi_imp2[,1])
      lt_significance <- colMeans(bt_condi_imp2 < bt_condi_imp2[,1])
      
    }
    
    conditional2 <- list(importance = condi_imp2, bt_imp = bt_condi_imp2,
                         benchmark = benchmark2, significance = significance, lt_significance = lt_significance)
    
    st2e_obj2$imp_scores <- list(joint = joint2, conditional = conditional2)
    
    list(st2e_obj, st2e_obj2)
    
  }else{
    st2e_obj
    
  }
  
}


# Print out the st2e selection
print.st2e<- function(st2e_obj, imp_type = c('conditional','joint'), sigthresh = 1-exp(-5)){
  cat('----------------------------------------------------\n')
  cat("Summary of the objective scores of selections in the ensemble:\n")
  print(summary(st2e_obj$obj_scores))
  cat('----------------------------------------------------\n')
  
  # Check the type of importance scores (and their significance) used to make statistically significant selections
  imp_type <- match.arg(imp_type)
  
  # Check whether importance scores were calculated
  if(is.null(st2e_obj$imp_scores)){
    cat("Importance scores were not found in the st2e object. Run calc_imp before to make statistically significant selections.")
    selected <- NULL
    
  }else{
    # Extract importance scores
    imp_obj <- st2e_obj$imp_scores[[imp_type]]
    
    imp <- imp_obj$importance[-1]
    sig <- imp_obj$significance[-1]
    
    cat('Selections are made based on', imp_type, 'importance scores \n')
    if(sum(sig>sigthresh)>0){
      cat('Predictors selected significantly more often than the pseudo-predictor are: \n')
      selected <- names(sig)[sig>sigthresh]
      cat(selected,'\n\n')
    }else{
      cat('Predictors with top 5 importancescores are: \n')
      selected <- names(imp)[rank(-imp, ties.method = 'min')<=5]
      cat(selected,'\n\n')
    }
    
  }
  
  invisible(list(selected = selected))
  
}


# Print out the stablemate selection
print.stablemate <- function(stbm_obj, pred_imp_type = c('joint','conditional'), 
                             stab_imp_type = c('conditional','joint'),  sigthresh = 1-exp(-5)){
  
  # Print out predictors selected as predictive
  cat('----------------------------------------------------\n')
  cat("Summary of the objective scores of the selections in the prediction ensemble:\n")
  print(summary(stbm_obj$prediction_ensemble$obj_scores))
  cat('----------------------------------------------------\n')
  pred_selected <- NULL
  
  # Check whether importance scores of predictivity selections were calculated
  if(is.null(stbm_obj$prediction_ensemble$imp_scores)){
    cat("Importance scores were not calculated on the prediction ensemble. Run calc_imp before to make statistically significant selections.")
    
  }else{
    # Check the type of importance scores (and their significance) used to make statistically significant selections
    pred_imp_type <- match.arg(pred_imp_type)
    
    # Extract the importance scores (and their significance) of predictive predictors
    pred_imp_obj <- stbm_obj$prediction_ensemble$imp_scores[[pred_imp_type]]
    pred_sel_imp <- pred_imp_obj$importance[-1]
    pred_sig <- pred_imp_obj$significance[-1]
    
    # Print out predictors selected as predictive
    cat('Selections of predictive predictors are made based on', pred_imp_type, 'importance scores \n')
    if(sum(pred_sig>sigthresh)>0){
      cat('Predictors selected significantly more often than the pseudo-predictor in the predictivity selection are: \n')
      pred_selected <- names(pred_sig)[pred_sig>sigthresh]
      cat(pred_selected,'\n\n')
    }else{
      cat('No predictor was found selected significantly more often than the pseudo-predictor in the predictivity selection. No selection will be made on stable or unstable predictors. \n\n')
    }
    
  }
  
  # Print out predictors selected as stable and predictive
  cat('----------------------------------------------------\n')
  cat("Summary of the objective scores of the selections in the stable ensemble:\n")
  print(summary(stbm_obj$stable_ensemble$obj_scores))
  cat('----------------------------------------------------\n')
  stab_selected = NULL;  nstab_selected = NULL
  
  # Check whether importance scores of stability selections were calculated
  if(is.null(stbm_obj$stable_ensemble$imp_scores)){
    cat("Importance scores were not calculated on the stable ensemble. Run calc_imp before to make statistically significant selections.")
    
  }else if(is.null(stbm_obj$prediction_ensemble$imp_scores)){
    cat("Selections on predictive and stable predictors cannot be made without calculating importance scores on the prediction ensemble.")
    
  }else{
    # Check the type of importance scores (and their significance) used to make statistically significant selections
    stab_imp_type <- match.arg(stab_imp_type)
    
    # Extract the importance scores (and their significance) of stable predictors
    stab_imp_obj <- stbm_obj$stable_ensemble$imp_scores[[stab_imp_type]]
    stab_sel_imp <- stab_imp_obj$importance[-1]
    stab_sig <- stab_imp_obj$significance[-1]
    nstab_sig <- stab_imp_obj$lt_significance[-1]
    
    # Print out predictors selected as predictive and stable
    if(!is.null(pred_selected)){
      cat('Selections of stale predictors are made based on', stab_imp_type, 'importance scores \n')
      if(sum(stab_sig>sigthresh)>0){
        stab_selected <- names(stab_sig)[stab_sig>sigthresh]
        stab_selected <- intersect(pred_selected, stab_selected)
        if(length(stab_selected) > 0){
          cat('Predictors selected as predictive and selected significantly more often than the pseudo-predictor in the stability selection are: \n')
          cat(stab_selected,'\n\n')
        }else{
          cat('No predictive predictor was found selected significantly more often than the pseudo-predictor in the stability selection. \n\n')
        }
      }else{
        cat('No predictor was found selected significantly more often than the pseudo-predictor in the stability selection. \n\n')
      }
      
      # Print out predictors selected as predictive but unstable
      if(sum(nstab_sig>sigthresh)>0){
        nstab_selected <- names(nstab_sig)[nstab_sig>sigthresh]
        nstab_selected <- intersect(pred_selected, nstab_selected)
        if(length(nstab_selected) > 0){
          cat('Predictors selected as predictive and selected significantly less often than the pseudo-predictor in the stability selection are: \n')
          cat(nstab_selected,'\n\n')
        }else{
          cat('No predictive predictor was found selected significantly less often than the pseudo-predictor in the stability selection. \n\n')
        } 
      }else{
        cat('No predictor was found selected significantly less often than the pseudo-predictor in the stability selection. \n\n')
      }
      
    }
    
  }
  
  invisible(list(pred_selected = pred_selected, stab_selected = stab_selected,  nstab_selected = nstab_selected))
  
}


# Plot out the st2e selection
plot.st2e <- function(st2e_obj, imp_type = c('conditional','joint'), 
                      fill_by = NULL, colors = NULL, sigthresh = 1-exp(-5),
                      plot_density = FALSE, bt_prop = 0.1, base_size = 15,
                      show_labels = TRUE, labels = NULL, label_size = 10, box_size = 0.1){
  
  # Check the type of importance scores (and their significance) used to make statistically significant selections
  imp_type <- match.arg(imp_type)
  
  # Check whether importance scores were calculated
  if(is.null(st2e_obj$imp_scores)){
    cat("Importance scores were not found in the st2e object. Run calc_imp before to make statistically significant selections.")
  }
  
  # Extract importance scores
  imp_obj <- st2e_obj$imp_scores[[imp_type]]
  
  # Name of predictors
  pred_names <- names(imp_obj$importance)[-1]
  
  # The number of predictors
  P <- length(pred_names)
  
  # The number of bootstraps
  B <- nrow(imp_obj$bt_imp)
  
  
  ############### Code for label tuning ############### 
  # Whether to show column labels (x-ticks)
  if(show_labels){
    # User provided x label, It should be a vector named by the predictor names 
    # to be shown
    if(!is.null(labels)){
      if(is.null(names(labels))) names(labels) <- labels
      tmp <- rep('',P)
      names(tmp) <- pred_names
      tmp[names(labels)] <- labels
      labels <- tmp
      
      # By default, show the labels of selected predictors
    }else{
      tmp <- rep('',P)
      names(tmp) <- pred_names
      tmp[imp_obj$significance[-1]>sigthresh] <- pred_names[imp_obj$significance[-1]>sigthresh]
      labels <- tmp
      
    }
    
  }
  
  
  ############### Code for the boxplot ############### 
  # Plot the bootstrapped distribution of the differences in importance scores between the pseudo predictor and each true predictor
  if(plot_density){
    # Randomly choose p percent of boostrap samples to show
    if(!is.null(bt_prop)) bt_id <- sample(1:B, bt_prop * B) else bt_id <- 1:B
    
    bt_imp_dff <- (imp_obj$bt_imp[bt_id,-1] - imp_obj$bt_imp[bt_id,1]) %>% as.matrix() %>% as.data.frame() %>%
      reshape2::melt(value.name = 'imp_diff')
    
    # Make a boxplot of the differences in importance scores between the pseudo predictor and each true predictor
    g <- ggplot2::ggplot(data = bt_imp_dff) + ggplot2::ylim(-1.001,1.001)  + 
      ggplot2::theme_grey(base_size = base_size) + 
      ggplot2::geom_boxplot(ggplot2::aes(variable, imp_diff, color = fill_by, fill = fill_by), size = box_size) +
      ggplot2::xlab('Predictors') + ggplot2::ylab(paste('Differences in', imp_type,'importance scores (Pseudo vs True)')) + 
      ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')
    
    
    ############### Code for color tuning ############### 
    # Create default filling according to the final selection result.
    if(is.null(fill_by)){
      fill_by <- rep('Not Selected', P)
      fill_by[imp_obj$significance[-1]>sigthresh] <- 'Selected'
      names(fill_by) <- pred_names
      fill_by <- fill_by[as.character(bt_imp_dff$variable)]
      
      g <- g + ggplot2::scale_fill_manual('Selection',
                                          values = c('lightgrey','red'))+
        ggplot2::scale_color_manual('Selection', values = c('lightgrey','red'))
      
    }else{
      # User provided filling, should be a data frame of the same row length as
      # the column size of the ensemble (exclude the pseudo-predictor),
      # Rows should be labelled by predictor names. The column name will be 
      # used as the legend title
      legend_title <- colnames(fill_by)
      fill_by <- fill_by[as.character(bt_imp_dff$variable),1]
      
      # Color schemes
      if(is.null(colors)){
        g <- g + ggplot2::scale_fill_discrete(legend_title)+
          ggplot2::scale_color_discrete(legend_title)
        
      }else{
        g <- g + ggplot2::scale_fill_manual(legend_title,
                                            values = colors)+
          ggplot2::scale_color_manual(legend_title, values = colors)
        
      }
      
    }
    
    
    ############### Code for the barplot ############### 
  }else{
    g <- ggplot2::ggplot() + ggplot2::ylim(0,1) + ggplot2::theme_classic(base_size = base_size) +
      ggplot2::geom_col(ggplot2::aes(x = pred_names, y = imp_obj$importance[-1], fill = fill_by), size = box_size)+
      ggplot2::geom_hline(yintercept = imp_obj$importance[1], linetype = 'dashed') +
      ggplot2::xlab('Predictors') + 
      if(imp_type == 'conditional') ggplot2::ylab('Conditional mportance scores') else ggplot2::ylab('Joint importance scores')
    
    
    ############### Code for color tuning ############### 
    # Create default filling according to the final selection result.
    if(is.null(fill_by)){
      fill_by <- rep('Not Selected', P)
      fill_by[imp_obj$significance[-1]>sigthresh] <- 'Selected'
      
      g <- g + ggplot2::scale_fill_manual('Selection',
                                          values = c('lightgrey','red'))+
        ggplot2::scale_color_manual('Selection', values = c('lightgrey','red'))
      
    }else{
      # User provided filling, should be a data frame of the same row length as
      # the column size of the ensemble (exclude the pseudo-predictor),
      # Rows should be labelled by predictor names. The column name will be 
      # used as the legend title
      legend_title <- colnames(fill_by)
      fill_by <- fill_by[pred_names,1]
      
      # Color schemes
      if(is.null(colors)){
        g <- g + ggplot2::scale_fill_discrete(legend_title)+
          ggplot2::scale_color_discrete(legend_title)
        
      }else{
        g <- g + ggplot2::scale_fill_manual(legend_title,
                                            values = colors)+
          ggplot2::scale_color_manual(legend_title, values = colors)
        
      }
      
    }
    
  }
  
  g + ggplot2::theme(axis.ticks.x =  ggplot2::element_blank(), 
                     axis.text.x= ggplot2::element_text(size = label_size)) +
    ggplot2::scale_x_discrete(labels = labels)
  
}


# Plot out the stablemate selection
plot.stablemate <- function(stbm_obj, pred_imp_type = c('joint','conditional'), 
                            stab_imp_type = c('conditional','joint'),
                            sigthresh = 1-exp(-5),
                            color_by = NULL, colors = NULL,
                            base_size = 15, point_size  = 2,
                            show_labels = TRUE, labels = NULL, label_size = 3, parse = FALSE,
                            pred_cutoff = NULL, show_density_of = NULL){
  
  
  ############### Initialization ############### 
  
  # Check the type of importance scores (and their significance) used to make statistically significant selections
  pred_imp_type <- match.arg(pred_imp_type)
  stab_imp_type <- match.arg(stab_imp_type)
  
  # Check whether importance scores of predictivity selections were calculated
  if(is.null(stbm_obj$prediction_ensemble$imp_scores)){
    stop("Importance scores were not calculated on the prediction ensemble. Run calc_imp before to make statistically significant selections.")
  }
  
  # Check whether importance scores of stability selections were calculated
  if(is.null(stbm_obj$stable_ensemble$imp_scores)){
    stop("Importance scores were not calculated on the stable ensemble. Run calc_imp before to make statistically significant selections.")
  }
  
  
  ############### Calculate plotting elements ############### 
  # Extract the importance scores (and their significance) of predictive predictors
  pred_imp_obj <- stbm_obj$prediction_ensemble$imp_scores[[pred_imp_type]]
  pred_imp <- pred_imp_obj$importance[-1]
  pred_sig <- pred_imp_obj$significance[-1]
  
  # Extract the importance scores (and their significance) of stable predictors
  stab_imp_obj <- stbm_obj$stable_ensemble$imp_scores[[stab_imp_type]]
  stab_imp <- stab_imp_obj$importance[-1]
  stab_sig <- stbm_obj$stable_ensemble$imp_scores[['conditional']]$significance[-1]
  nstab_sig <- stbm_obj$stable_ensemble$imp_scores[['conditional']]$lt_significance[-1]
  
  # Log transform the significance of the importance scores in stability selection
  # which can be used to color predictors in the variable selection plot. 
  # Its a measure of how significant, either stable or unstable, each predictor is
  log_stab_sig <- 1 - (pmax(stab_sig, nstab_sig) - pmin(stab_sig, nstab_sig))
  log_stab_sig <- -log10(log_stab_sig)
  # Unstable predictors receive negative log-significance, whereas stable predictors receive
  # positive log-significance. The larger the absolute value is, the more stable (or unstable)
  # a predictors
  log_stab_sig[nstab_sig > stab_sig] <- -log_stab_sig[nstab_sig > stab_sig]
  log_stab_sig[nstab_sig == stab_sig] <- 0
  
  # Make the final selection based on the significance threshold provided
  is_pred <- pred_sig > sigthresh # Whether a predictor is predictive
  is_stab <- rep('Non-significant',length(stab_sig))
  is_stab[stab_sig > sigthresh] <- 'Stable'
  is_stab[nstab_sig > sigthresh] <- 'Unstable'
  is_stab <- factor(is_stab, levels = c('Stable','Non-significant','Unstable'))
  
  # Set the upper bound of log-pvalue
  B <- nrow(stab_imp_obj$bt_imp) # The size of the ensemble
  p_upbd <- log10(B+1)
  log_stab_sig[log_stab_sig == Inf] <- p_upbd
  log_stab_sig[log_stab_sig == -Inf] <- -p_upbd
  
  # Set up a color range, where the predictors with significance higher than the 
  # threshold will be marked with specific color coding for stability
  q <- (p_upbd - min(-log10(sigthresh), p_upbd))/(2*p_upbd)
  
  # Set up predictors to color
  pred_names <- names(pred_imp)
  alpha <- rep(0, length(pred_names))
  names(alpha) <- pred_names
  if(!is.null(pred_cutoff)){
    alpha[pred_imp>pred_cutoff] <- 1
  }else{
    alpha[pred_sig>sigthresh] <- 1
    
  }
  
  # Set up labels
  if(show_labels){
    # User provided x label, It should be a vector named by the predictor names 
    # to be shown
    if(!is.null(labels)){
      label_which <- ifelse(is.null(names(labels)), labels, names(labels)) 
      
      # By default, show the labels of predictive predictors
    }else if(is.null(pred_cutoff)){
      label_which <- labels <- pred_names[pred_sig>sigthresh]
      
    }else{
      label_which <- labels <- pred_names[pred_imp>pred_cutoff]
      
    }
    
  }
  
  
  ############### Make plot ############### 
  g <- ggplot2::ggplot()
  
  # Draw bootstrap quantiles of the pseudo-predictor
  if(stab_imp_type == 'conditional'){
    xintercept <- quantile(pred_imp_obj$bt_imp[,1], c(1-sigthresh, sigthresh))
    yintercept <- quantile(stab_imp_obj$bt_imp[,1], c(1-sigthresh, sigthresh))
    
    g <- g + ggplot2::geom_vline(xintercept = xintercept, linetype = 'dashed', col = 'blue', alpha = 0.5) + 
      ggplot2::geom_hline(yintercept = yintercept, linetype = 'dashed', col = 'blue', alpha = 0.5)+
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= 0, ymax = 1), fill = 'green', alpha = 0.2) +
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= yintercept[2], ymax = 1), fill = 'green', alpha = 0.2) +
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= 0, ymax = yintercept[1]), fill = 'blue', alpha = 0.1)
    
  }else{
    xintercept <- quantile(stbm_obj$prediction_ensemble$imp_scores$conditional$bt_imp[,1], c(1-sigthresh, sigthresh))
    slope <- quantile(stbm_obj$stable_ensemble$imp_scores$conditional$bt_imp[,1], c(1-sigthresh, sigthresh))
    yintercept <- c(0,0)
    
    g <- g + ggplot2::geom_vline(xintercept = xintercept, linetype = 'dashed', col = 'blue', alpha = 0.5) +
      ggplot2::geom_abline(slope = slope, intercept = yintercept, linetype = 'dashed', col = 'blue', alpha = 0.5)+
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= 0, ymax = 1), fill = 'green', alpha = 0.2) +
      ggplot2::geom_ribbon(ggplot2::aes(x = seq(xintercept[2], 1, length.out = 1000), 
                                        ymin = seq(xintercept[2], 1, length.out = 1000) * slope[2] + yintercept[2],
                                        ymax = 1), fill = 'green', alpha = 0.2)+
      ggplot2::geom_ribbon(ggplot2::aes(x = seq(xintercept[2], 1, length.out = 1000), 
                                        ymax = seq(xintercept[2], 1, length.out = 1000) * slope[1] + yintercept[1],
                                        ymin = 0), fill = 'blue', alpha = 0.1)
  }
  
  
  ############### Density plot ############### 
  # Plot the 2d density of the bootstrapped importance scores of user selected predictors and the pseudo-predictor
  if(!is.null(show_density_of)){
    # Extract the bootstrapped importance scores of user selected predictors 
    message(cat('Chosen to plot the bootstrapped selection probabilities of', show_density_of))
    bt_pred_imp <- pred_imp_obj$bt_imp[,show_density_of,drop=FALSE]
    bt_stab_imp <- stab_imp_obj$bt_imp[,show_density_of,drop=FALSE]
    
    getLevel <- function(x,y,prob= 2*sigthresh-1) {
      kk <- MASS::kde2d(x,y)
      dx <- diff(kk$x[1:2])
      dy <- diff(kk$y[1:2])
      sz <- sort(kk$z)
      c1 <- cumsum(sz) * dx * dy
      approx(c1, sz, xout = 1 - prob)$y
    }
    
    # Draw the density line of the pseudo-predictor
    L95 <- getLevel(pred_imp_obj$bt_imp[,1], stab_imp_obj$bt_imp[,1])
    g <- g +  ggplot2::geom_density_2d(ggplot2::aes(pred_imp_obj$bt_imp[,1], stab_imp_obj$bt_imp[,1]), color = 'black', breaks = L95)
    
    for(i in show_density_of){
      # Draw the density line of the ith user selected predictor
      L95 <- getLevel(bt_pred_imp[,i],bt_stab_imp[,i])
      g <- g +  ggplot2::geom_density_2d(ggplot2::aes_string(bt_pred_imp[,i], bt_stab_imp[,i]), color = 'purple', breaks = L95)
      
      # If significance was evaluated at per-predictor level (i.e, each true predictor is 
      # benchmarked against the bootstrap distribution of the pseudo-predictor selected 
      # under the pool of the true predictor)
      if(!is.null(pred_imp_obj$benchmark) | !is.null(stab_imp_obj$benchmark)){
        if(is.null(pred_imp_obj$benchmark)){
          benchmark_x <- pred_imp_obj$bt_imp[,1]
          
        }else{
          benchmark_x <- pred_imp_obj$benchmark[,i]
          
        }
        
        if(is.null(stab_prob_obj$benchmark)){
          benchmark_y <- stab_imp_obj$bt_imp[,1]
          
        }else{
          benchmark_y <- stab_imp_obj$benchmark[,i]
          
        } 
        
        # Draw the density line of the pseudo-predictor selected under the pool of the true predictor
        L95 <- getLevel(benchmark_x, benchmark_y)
        g <- g +  ggplot2::geom_density_2d(ggplot2::aes_string(benchmark_x, benchmark_y), color = 'black', breaks = L95)
        g <- g +  ggplot2::annotate('text', x = mean(benchmark_x), y = mean(benchmark_y), label = i,
                                    color = 'black',size = label_size)
        
      }
      
    }
    
  } 
  
  
  ############### The main plot ############### 
  g <- g + ggplot2::scale_x_continuous(limits = c(-0.001,1.001), expand = c(0, 0.01)) +
    ggplot2::scale_y_continuous(limits = c(-0.001,1.001), expand = c(0, 0.01)) + ggplot2::theme_grey(base_size = base_size) +
    ggplot2::geom_point(ggplot2::aes(pred_imp, stab_imp, shape = is_stab),
                        color = 'lightgrey',fill = 'lightgrey', alpha = 1-alpha, size = point_size)+
    ggrepel::geom_text_repel(ggplot2::aes(x = pred_imp[label_which], y = stab_imp[label_which], label = labels), 
                             segment.color = 'grey', force = 1,max.overlaps = 20, size = label_size, segment.linetype = 'dashed', 
                             parse = parse) + 
    ggplot2::scale_shape_manual('Selection', values = c('Stable' = 24, 'Non-significant' = 21, 'Unstable' = 25))
  
  g <- g + ggplot2::geom_point(ggplot2::aes(pred_imp, stab_imp,
                                            fill = color_by, shape = is_stab, col = color_by), alpha = alpha, size = point_size) 
  
  
  if(is.null(color_by)){
    color_by <- is_stab
    g <- g + ggplot2::scale_fill_manual('Selection', values = c('Stable' = "blue", 'Non-significant' = "orange", 'Unstable' = "red")) + 
      ggplot2::scale_color_manual('Selection', values = c('Stable' = "blue", 'Non-significant' = "orange", 'Unstable' = "red"))
    
  }else if(length(color_by) == 1 & color_by[1] == 'p-val'){
    color_by <- log_stab_sig
    g <- g + ggplot2::scale_fill_gradientn('Signed Log10 P-val',colours =  c("red", "orange","orange","blue"),
                                           values = c(0, q, 1-q, 1), limits = c(-p_upbd,p_upbd))+
      ggplot2::scale_color_gradientn('Signed Log10 P-val',colours =  c("red", "orange", "orange","blue"),
                                     values = c(0, q, 1-q, 1), limits = c(-p_upbd,p_upbd))
    
  }else{
    # User provided coloring, should be a data frame of the same row length as
    # the column size of the ensemble (exclude the pseudo-predictor),
    # Rows should be labelled by predictor names. The column name will be 
    # used as the legend title
    legend_title <- colnames(color_by)
    color_by <- color_by[as.character(pred_names),1]
    
    # Color schemes
    if(is.null(colors)){
      g <- g + ggplot2::scale_fill_discrete(legend_title)+
        ggplot2::scale_color_discrete(legend_title)
      
    }else{
      g <- g + ggplot2::scale_fill_manual(legend_title,
                                          values = colors)+
        ggplot2::scale_color_manual(legend_title, values = colors)
      
    }
    
  }
  
  x <- "Importance: Predictive"
  if(stab_imp_type == 'conditional'){
    y <- "Importance: Stable"
  }else{
    y <- "Importance: Stable and Predictive"
  }
  
  g + ggplot2::xlab(x) + ggplot2::ylab(y)
  
}


# Lasso pre-filtering
lasso_prefilt <- function(Y, X, p, env = NULL, K = 100, ncore = NULL, drop_pred = TRUE,
                          par_method = c('SNOW','MC','MPI'), verbose = TRUE,
                          chunk_size = ceiling(K/ncore), group_by = c('env', 'Yenv', 'none'), ...){
  # If environment variable is not provided, create a pseudo environment that includes all samples
  if(is.null(env)) env <- rep('pseudo_environment', N)
  # The number of samples
  N <- nrow(X)
  # The number of predictors
  P <- ncol(X)
  # Predictor names
  pred_names <- colnames(X)
  
  # Define sample group
  group_by <- match.arg(group_by)
  if(group_by == 'env'){
    group <- env
    
    # Used for logistic regression
  }else if(group_by == 'Yenv'){
    if(length(unique(Y)) != 2) warning('Y has more than two unique values, group samples by Yenv may not be valid')
    group <- factor(paste(Y, env, sep = '_'))
    
  }else{
    group <- 1
    
  }
  
  # Define one iteration of lasso selection
  one_lasso <- function(){
    # Sample data in half
    idx <- sample(1:N, ceiling(N/2)) 
    # The number of samples used to run lasso selection
    n <- length(idx)
    
    # Define sample weights
    if(group_by != 'none'){
      w <- 1/(table(group[idx])[group[idx]]) 
      w <- as.numeric(w/sum(w)*n)
      
      # Run Lasso on subsampled data
      suppressWarnings(mod <- glmnet::glmnet(X[idx,], Y[idx], weights = w, ...))
    }else{
      suppressWarnings(mod <- glmnet::glmnet(X[idx,], Y[idx], ...))
      
    }
    
    # Indicator matrix that record the lasso selection path
    path_mat <- (mod$beta != 0) 
    # Calculate at which step each predictor is first selected    
    first_entrance <- apply(path_mat, 1, which.max) 
    # If not selected in the entire path, a predictor is set to be selected at a order equals to infinity
    first_entrance[which(apply(path_mat, 1, sum) == 0)] <- Inf
    # Calculate the number of predictors that are not selected by the entire selection path    
    num_noinf <- sum(first_entrance != Inf) 
    # Select the first p predictors entering the path
    sel <- first_entrance <= sort(first_entrance)[min(p,num_noinf)]
  }
  
  # Run lasso selection
  if(is.null(ncore)){
    lasso_sel <- matrix(0, ncol = ncol(X), nrow = K)
    colnames(lasso_sel) <- colnames(X)
    
    # Retrieved from @https://www.dummies.com/article/technology/programming-web-design/r/how-to-generate-your-own-error-messages-in-r-175112/
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = K,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
    
    for(i in 1:K){
      pb$tick() # Timing progress bar
      lasso_sel[i,] <- one_lasso()
    }
    
    # Run lasso with doParallel, which should work on all platform
  }else{
    par_method <- match.arg(par_method)
    `%dopar%` <- getFromNamespace("%dopar%", "foreach")
    
    # Set up cluster
    ############### Run st2e with doParallel, which should work on all platforms ###############
    # SNOW is preferred on a Windows Machine
    if(par_method == 'SNOW'){
      cl<- parallel::makeCluster(ncore)
      doParallel::registerDoParallel(cl)
      
      
    }else if(par_method == 'MC'){
      ############### Run st2e with MC ###############
      # On a non-Windows Machine (e.g. Mac OS, Linux), MC essentially performs the multicore
      # functionality in the Parallel package
      doParallel::registerDoParallel(cores = ncore)
      
      
    }else{
      ############### Run st2e with MPI ############### 
      # MPI should be used when running on a HPC cluster, and
      # should be ran together with command line mpirun in shell scripts.
      if(is.null(doMPI::getDoMpiCluster())){
        cl <- doMPI::startMPIcluster(ncore)
        doMPI::registerDoMPI(cl)
        
      }
      
    }
    
    # Foe each loop
    if(chunk_size == 1){
      res <- foreach::foreach(i=1:K, .verbose = verbose) %dopar% {
        one_lasso()      
        
      }
      
    }else{
      idx_chunk <- split(1:K, ceiling(seq_along(1:K)/chunk_size))
      
      res <- foreach::foreach(i=1:length(idx_chunk), .verbose = verbose) %dopar% {
        sel <- matrix(nrow = length(idx_chunk[[i]]), ncol = ncol(X))
        for(j in 1:length(idx_chunk[[i]])){
          sel[j,] <- one_lasso()
          
        }
        sel
        
      }
      
    }
    
    lasso_sel <- do.call(rbind, res)
    colnames(lasso_sel) <- colnames(X)
    if(par_method == 'SNOW') parallel::stopCluster(cl)
    
  }
  
  if(drop_pred) lasso_sel <- lasso_sel[,colSums(lasso_sel)>0]
  lasso_sel <- cbind(pseudo_predictor = 1, lasso_sel)
  lasso_sel <- Matrix::Matrix(lasso_sel)
  attributes(lasso_sel)$stbm_lasso <- "TRUE"
  lasso_sel
  
}


# Given objective scores of multiple proposals, defines how to select a proposal
sel_min <- function(scores){
  which.min(scores)
}


# Make prediction using st2e
predict.st2e <- function(st2e_obj, X, prune = TRUE, prefilt_scores = NULL, pred_fun = pred_ols){
  K <- nrow(st2e_obj$ensemble)
  # Calculate the weights of st2e_obj selections if prune the ensemble
  if(prune) q_sel <- rank(-st2e_obj$obj_scores, ties.method = 'average')/K else q_sel <- rep(1,K)
  # Calculate the weights of prior selections if prune the ensemble 
  if(prune & !is.null(prefilt_scores)) q_prefilt <- rank(-prefilt_scores, ties.method = 'average')/K else q_prefilt <-  rep(1,K)
  q_sel <- q_sel * q_prefilt
  
  # Build prediction ensemble
  prediction_ensemble <- matrix(0, ncol = nrow(X), nrow = K)
  pred_selected <- names(which(colSums(st2e_obj$ensemble[,-1])>0)) # Check whether all predictors that were used to build models are also in the new data
  if( !all( pred_selected %in% colnames(X) )) stop('At least one predictor in the ensemble was not found in the data')
  
  #Each model makes prediction
  for(i in 1:K){
    pred_selected <- names(which(st2e_obj$ensemble[i,-1]==1))
    prediction_ensemble[i,] <- pred_fun(st2e_obj$model[[i]], X[, pred_selected, drop =F])
  }
  
  aggregated_prediction <- colSums(prediction_ensemble*q_sel)/sum(q_sel) 
  
  attr(aggregated_prediction, which = 'ensemble') <- prediction_ensemble
  
  aggregated_prediction
}


# Make prediction using stablemate
predict.stablemate <- function(stbm_obj, X, assay = c('stable_ensemble','prediction_ensemble'),
                               prune = TRUE, prefilt_scores = NULL, pred_fun = pred_ols){
  assay <- match.arg(assay)
  
  # Check whether the stability selection is based on the predictivity selection
  if(stbm_obj$matched){
    K <- nrow(stbm_obj$stable_ensemble$ensemble)
    
    if(prune){
      # Calculate the weights of stability selections if prune the ensemble (In stablemate, its the ranking quantiles of pss scores)
      q_sel <- rank(-stbm_obj$stable_ensemble$obj_scores, ties.method = 'average')/K 
      # Calculate the weights of st2e_obj2 selections if prune the ensemble (In stablemate, its the ranking quantiles of bic scores)
      q_sel2 <- rank(-stbm_obj$prediction_ensemble$obj_scores, ties.method = 'average')/K 
      # Calculate the weights of prior selections if prune the ensemble (In stablemate, its the ranking quantiles of the objective scores of the pre-filtering method applied)
      if(!is.null(prefilt_scores)) q_prefilt <- rank(-prefilt_scores, ties.method = 'average')/K else q_prefilt <-  rep(1,K)
      
    }else{
      q_sel <- q_sel2 <- q_prefilt <- rep(1,K)
      
    }
    
    # Weight selections by whether they are both important for the pre-filtering and predictivity selection
    q_sel2 <- q_sel2*q_prefilt
    # Weight selections by whether they are both important for the pre-filtering and predictivity selection and stability selection
    q_sel <- q_sel*q_sel2
    
    # Build prediction ensemble
    ensemble <- matrix(0, ncol = nrow(X), nrow = K)
    pred_selected <- names(which(colSums(stbm_obj[[assay]]$ensemble[,-1])>0)) # Check whether all predictors that were used to build models are also in the new data
    if( !all( pred_selected %in% colnames(X) )) stop('At least one predictor in the ensemble was not found in the data')
    
    #Each model makes prediction
    for(i in 1:K){
      pred_selected <- names(which(stbm_obj[[assay]]$ensemble[i,-1]==1))
      ensemble[i,] <- pred_fun(stbm_obj[[assay]]$model[[i]], X[, pred_selected, drop =F])
    }
    
    aggregated_prediction <- if(assay == 'stable_ensemble') colSums(ensemble*q_sel)/sum(q_sel) else
      colSums(ensemble*q_sel2)/sum(q_sel2)
    
    attr(aggregated_prediction, which = 'ensemble') <- ensemble
    aggregated_prediction
    
  }else{
    predict.st2e(stbm_obj[[assay]], X, prune = prune, prefilt_scores = prefilt_scores, pred_fun = pred_fun)
    
  }
  
}


# Default regression and objective functions
obj_bic <- function(Y, X, env, pred_in){
  n_sample <- length(Y)
  env_tab <- table(env)
  n_in <- sum(pred_in)
  
  # Calculate sample weights
  if(length(env_tab) == 1){
    weight <- FALSE
    w <- rep(1, n_sample)
  }else{
    weight <- TRUE
    w <- 1/(env_tab[env]) # Weight samples by the (inverse) sizes of their environment
    w <- as.numeric(w/sum(w)*n_sample)
  }
  
  if(n_in>0){
    mod <- if(weight) lm.wfit(cbind(1,X[,pred_in==1]), Y, w) else lm.fit(cbind(1,X[,pred_in==1]), Y)
    res <- mod$residuals
    p <- mod$rank + 1
  }else{
    mod <- mean(w*Y)
    res <- Y - mod
    p <- 2
  }
  
  ll <- .5* (sum(log(w)) - n_sample * (log(2 * pi) + 1 - log(n_sample) +
                                         log(sum(w*res^2))))
  
  p  * log(n_sample) - 2*ll
}
reg_ols <- function(Y, X, env, pred_in){
  n_sample <- length(Y)
  env_tab <- table(env)
  n_in <- sum(pred_in)
  
  # Calculate sample weights
  if(length(env_tab) == 1){
    weight <- FALSE
    w <- rep(1, n_sample)
  }else{
    weight <- TRUE
    w <- 1/(env_tab[env]) # Weight samples by the (inverse) sizes of their environment
    w <- as.numeric(w/sum(w)*n_sample)
  }
  
  if(n_in>0){
    mod <- if(weight) lm.wfit(cbind(1,X[,pred_in==1]), Y, w) else lm.fit(cbind(1,X[,pred_in==1]), Y)
    coef <- mod$coefficients
  }else{
    coef <- mean(w*Y)
  }
  
  out <- rep(0, ncol(X)+1)
  names(out) <- c('(Intercept)', colnames(X))
  out[colnames(X)[pred_in==1]] <- coef[-1]
  out[1] <- coef[1]
  out
}
obj_pss <- function(Y, X, env, pred_in){
  n_sample <- length(Y)
  env_tab <- table(env)
  n_in <- sum(pred_in)
  
  # Calculate sample weights
  if(length(env_tab) == 1){
    weight <- FALSE
    w <- rep(1, n_sample)
  }else{
    weight <- TRUE
    w <- 1/(env_tab[env]) # Weight samples by the (inverse) sizes of their environment
    w <- as.numeric(w/sum(w)*n_sample)
  }
  
  pss <- 0
  for(e in names(env_tab)){
    
    if(n_in>0){
      mod <- if(weight) lm.wfit(cbind(1,X[env != e, pred_in==1]), Y[env != e], w[env != e]) else
        lm.fit(cbind(1,X[env != e, pred_in==1]), Y[env != e])
      coef <- mod$coefficients
      coef[is.na(coef)] <- 0
      if (length(env_tab) > 2) {
        Yhat <- tcrossprod(X[env == e, pred_in==1], t(coef[-1])) + coef[1]
        ss <- (Yhat - Y[env == e])^2
        pss <- pss + sum(w[env == e] * ss)
        
        mod <- if(weight) lm.wfit(cbind(1,X[env == e, pred_in==1]), Y[env == e], w[env == e]) else
          lm.fit(cbind(1,X[env == e, pred_in==1]), Y[env == e])
        Yhat <- mod$fitted.values
        ss <- (Yhat - Y[env == e])^2
        pss <- pss + sum(w[env == e] * ss)
        
      }else{
        Yhat <- tcrossprod(X[, pred_in==1], t(coef[-1])) + coef[1]
        ss <- (Yhat - Y)^2
        pss <- pss + sum(w * ss)
        
      }
      
    }else{
      Yhat <- sum(w[env != e]*Y[env != e])/sum(w[env != e])
      
      if (length(env_tab) > 2) {
        ss <- (Yhat - Y[env == e])^2
        pss <- pss + sum(w[env == e] * ss)
        
        Yhat <-  sum(w[env == e]*Y[env == e])/sum(w[env == e])
        ss <- (Yhat - Y[env == e])^2
        pss <- pss + sum(w[env == e] * ss)
        
      }else{
        ss <- (Yhat - Y)^2
        pss <- pss + sum(w * ss)
        
      }
      
    }
    
  }
  pss
  
}
pred_ols <- function(model, X){
  model[1] + X %*% model[colnames(X)]
}


# Optional regression and objective functions
obj_bic_logit <- function(Y, X, env, pred_in){
  n_sample <- length(Y)
  n_in <- sum(pred_in)
  group <- factor(paste(Y, env, sep = '_'))
  group_tab <- table(group)
  w <- 1/(group_tab[group]) # Weight samples by the (inverse) sizes of their environment
  w <- as.numeric(w/sum(w)*n_sample)
  
  if(n_in>0){
    mod <- arm::bayesglm.fit(cbind(1,X[,pred_in==1]), Y, w,
                             family = binomial(),control = glm.control(maxit = 500)) 
    mu <- mod$fitted.values
    p <- mod$rank
    
  }else{
    mu <- mean(w*Y)
    p <- 1
    
  }
  ll <- sum(w*dbinom(Y, 1, mu, log =T))
  
  p  * log(n_sample) - 2*ll
}
reg_logit <- function(Y, X, env, pred_in){
  n_sample <- length(Y)
  n_in <- sum(pred_in)
  group <- factor(paste(Y, env, sep = '_'))
  group_tab <- table(group)
  w <- 1/(group_tab[group]) # Weight samples by the (inverse) sizes of their environment
  w <- as.numeric(w/sum(w)*n_sample)
  
  if(n_in>0){
    mod <- arm::bayesglm.fit(cbind(1,X[,pred_in==1]), Y, w,
                             family = binomial(),control = glm.control(maxit = 500)) 
    coef <- mod$coefficients
    
  }else{
    coef <- mean(w*Y)
    
  }
  
  out <- rep(0, ncol(X)+1)
  names(out) <- c('(Intercept)', colnames(X))
  out[colnames(X)[pred_in==1]] <- coef[-1]
  out[1] <- coef[1]
  out
}
obj_pss_logit <- function(Y, X, env, pred_in){
  n_sample <- length(Y)
  n_in <- sum(pred_in)
  env_tab <- table(env)
  group <- factor(paste(Y, env, sep = '_'))
  group_tab <- table(group)
  w <- 1/(group_tab[group]) # Weight samples by the (inverse) sizes of their environment
  w <- as.numeric(w/sum(w)*n_sample)
  ilogit <- binomial()$linkinv
  
  pss <- 0
  for(e in names(env_tab)){
    if(n_in>0){
      mod <- arm::bayesglm.fit(cbind(1,X[env != e, pred_in==1]), Y[env != e], w[env != e],
                               family = binomial(),control = glm.control(maxit = 500)) 
      coef <- mod$coefficients
      coef[is.na(coef)] <- 0
      if (length(env_tab) > 2) {
        Yhat <- ilogit(tcrossprod(X[env == e, pred_in==1], t(coef[-1])) + coef[1])
        
        ss <- -2 * dbinom(Y[env == e],1,Yhat,log = T)
        pss <- pss + sum(w[env == e] * ss)
        
        mod <- arm::bayesglm.fit(cbind(1,X[env == e, pred_in==1]), Y[env == e], w[env == e],
                                 family = binomial(),control = glm.control(maxit = 500)) 
        Yhat <- mod$fitted.values
        ss <- -2 * dbinom(Y[env == e],1,Yhat,log = T)
        pss <- pss + sum(w[env == e] * ss)
        
      }else{
        Yhat <- ilogit(tcrossprod(X[, pred_in==1], t(coef[-1])) + coef[1])
        ss <-  -2 * dbinom(Y,1,Yhat,log = T)
        pss <- pss + sum(w * ss)
        
      }
      
    }else{
      mod <- arm::bayesglm.fit(matrix(1,sum(env != e)), Y[env != e], w[env != e],
                               family = binomial(), control = glm.control(maxit = 500)) 
      Yhat <- ilogit(mod$coefficients)
      if (length(env_tab) > 2) {
        ss <- -2 * dbinom(Y[env == e],1,Yhat,log = T)
        pss <- pss + sum(w[env == e] * ss)
        
        mod <- arm::bayesglm.fit(matrix(1,sum(env == e)), Y[env == e], w[env == e],
                                 family = binomial(), control = glm.control(maxit = 500)) 
        Yhat <- ilogit(mod$coefficients)
        ss <- -2 * dbinom(Y[env == e],1,Yhat,log = T)
        pss <- pss + sum(w[env == e] * ss)
        
      }else{
        ss <- -2 * dbinom(Y,1,Yhat,log = T)
        pss <- pss + sum(w * ss)
        
      }
      
    }
  }
  pss
}
pred_logit <- function(model, X){
  ilogit <- function(x) exp(x)/(1+exp(x))
  ilogit(model[1] + X %*% model[colnames(X)])
}
