# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build pre-filtering ensemble by a random Lasso procedure.
#'
#' Build pre-filtering ensemble by a random Lasso procedure that is similar to
#' \href{https://doi.org/10.1111/j.1467-9868.2010.00740.x}{Stability selection, Meinshausen et al. (2010, The Journal of the Royal Statistical Society, Series B)})
#' For each st2 run, we randomly sample one half of samples to select the top \code{p} predictors with Lasso.
#' The pre-filtered predictors is then set as the predictor pool. The advantages are of two folds.
#' First, across the different re-sampling runs, the top \code{p} predictors are expected to differ thus
#' enabling to cover a large and diverse range of predictors in our overall search.
#' Second, we improve the stability of Lasso pre-screening when subjected to sample perturbation.
#'
#' @param Y A response vector or matrix depending on the objective function. It should be a vector if the function is used with
#' StableMate's default objective.
#' @param X A predictor matrix with rows representing samples and columns representing predictors. The columns must be named.
#' @param p Numerical; Pre-filter size for each random Lasso procedure.
#' @param K Numerical; Ensemble size.
#' @param drop_pred Logical; If true, remove predictors from X that are not in the predictor pool to save space and time.
#' @param ncore Numerical; Numerical; If greater than 0. Parallel computing is enabled.
#' @param par_method; Parallel computing method. SNOW is preferred on Windows local machines, MC is preferred on non-Windows local machines.
#' MPI is preferred on distributed computing cloud.
#' @param verbose Logical; If true, multi-core verbose will be printed.
#' @param chunk_size Numerical; The size of task chunks (How many st2 runs are performed one each task).
#' @param group_by Either 'env', 'Yenv' or 'none' to weight samples by the inverse sizes of their environments,
#' their Y class-environments combined if binomial regression, or no weighting respectively.
#' @param ... Arguments passed to \code{\link[glmnet]{glmnet}}.
#'
#' @export
#'
#' @rdname lasso_prefilt
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
