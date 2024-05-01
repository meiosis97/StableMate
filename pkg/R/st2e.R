# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Not exported satellite functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generator of beta-binomial distribution
rbeta_binom <- function(n, size, shape1, shape2) rbinom(n, size, prob = rbeta(n, shape1, shape2))


# Probability mass function of beta-binomial distribution
dbeta_binom <- function(x, size, shape1, shape2) choose(size,x) * beta(x + shape1, size - x + shape2) / beta(shape1, shape2)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Stochastic stepwise variable selection ensemble
#'
#' The main st2e function for stochastic stepwise variable selection proposed by
#' \href{https://www.tandfonline.com/doi/full/10.1080/10618600.2012.679223}{Lu Xin &Mu Zhu (2010, Journal of Computational and Graphical Statistics)}.
#'
#' @param Y A response vector or matrix depending on the objective function. It should be a vector if the function is used with
#' StableMate's default objective.
#' @param X A predictor matrix with rows representing samples and columns representing predictors. The columns must be named.
#' @param obj_fun The function that calculates objective scores for any given set of predictors. See our default objective function to learn the syntax. By
#' default the objective is to minimize linear regression BIC scores. See \code{\link[StableMate]{obj_fun}} for the other options and how to customize this function.
#' @param sel_fun The function that selects the best set of st2e proposal given multiple proposals and their objective scores. By
#' default we select the proposal with the minimum objective score. See \code{\link[StableMate]{sel_fun}} for how to customize this function.
#' @param reg_fun The regression function. By default the we perform ordinary least square regression.
#'  See \code{\link[StableMate]{reg_fun}} for the other options and how to customize this function.
#' @param pred_pool A vector or a matrix that restricts the variable selection space. If a vector is given, it must be a vector of characters that is a subset of
#' the column names of X. Section of "pseudo_predictor" can be enabled by adding a "pseudo_predictor" into \code{pred_pool}. When a matrix is given, it must be an indicator matrix
#' of \code{K} rows with column names being a subset of the column names of x. Section of "pseudo_predictor" can be enabled by adding a "pseudo_predictor" column into \code{pred_pool}.
#' @param start_set A vector or a matrix that controls the starting point of st2e selection. If a vector is given, it must be a vector of characters that is a subset of \code{pred_pool}
#' When a matrix is given, it must be an indicator matrix of \code{K} rows with column names being a subset of the column names of \code{pred_pool}. \code{start_set} must be smaller or
#' equal to \code{pred_pool} in every elements of their matching column. By default, we start the selection from a random predictor set under \code{pred_pool}.
#' @param max_size Numerical; Any step of selection will not consider a predictor set that is larger than \code{max_size}.
#' @param alpha Numerical; A Parameter of beta distribution for sampling proposal size.
#' @param beta Numerical; A Parameter of beta distribution for sampling proposal size.
#' @param lambda Numerical; A scaling factor that controls the size of binomial distribution for sampling proposal size.. Should be a value between 0 and 1.
#' If 1, then the proposal size can be as large as the size of the maximum possible proposal.
#' @param a Numerical; A scaling factor that controls the number of proposal to sample, which contains a term where \code{a} is multiplied by the size of the maximum possible proposal
#' to the power of \code{t}.
#' @param b Numerical; A scaling factor that controls the number of proposal to sample, which contains a term where \code{b} is multiplied by a quality that measures the number of
#' all possible unique proposal that can be made.
#' @param t Numerical; A scaling factor that controls the number of proposal to sample, which contains a term where the size of the maximum possible proposal is rasied to the power of
#' \code{t} and multiplied by \code{a}.
#' @param K Numerical; Ensemble size.
#' @param do_switch Logical; If true, switch steps that swap predictors that are in and out the model are allowed.
#' @param ret_mod Logical; If true, return the prediction ensemble
#' @param ret_imp Logical; If true, return the importance score.
#' @param calc_imp_ctrl A list that controls importance score calculation. See the document of \code{\link[StableMate]{calc_imp}} for details.
#' @param ncore Numerical; Numerical; If greater than 0. Parallel computing is enabled.
#' @param par_method; Parallel computing method. SNOW is preferred on Windows local machines, MC is preferred on non-Windows local machines.
#' MPI is preferred on distributed computing cloud.
#' @param verbose Logical; If true, multi-core verbose will be printed.
#' @param chunk_size Numerical; The size of task chunks (How many st2 runs are performed one each task).
#' @param fun_export A vector of names of functions to be imported by \code{\link[foreach]{foreach}} environment.
#' @param pkg_export A vector of names of packages to be imported by \code{\link[foreach]{foreach}} environment.
#' @param drop_pred Logical; If true, remove predictors from X that are not in the predictor pool to save space and time.
#'
#' @return A \code{st2e} object.
#'
#' @export
#' @name st2e
#' @rdname st2e
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
    if(class(pred_pool) == "lasso_ensemble" & is.null(calc_imp_ctrl$scale)){
      calc_imp_ctrl$scale <- TRUE

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
      if(!"doMPI" %in% installed.packages()[,1]){
        stop("Package doMPI not detected. Please install from CRAN before running st2e.")
      }

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
    default <- list(prune = T, prefilt_scores = NULL, scale = F, pooled = T, B = 5000)

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


#' The function for making forward selections in ST2
#'
#' @export
#' @rdname st2e
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


#' The function for making backward selections in ST2
#'
#' @export
#' @rdname st2e
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


#' The function for making switching of predictors in ST2
#'
#' @export
#' @rdname st2e
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
