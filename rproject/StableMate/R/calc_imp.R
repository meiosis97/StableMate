# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Calculate importance scores on a \code{st2e} object.
#'
#' \code{\link[StableMate]{st2e}} returns a \code{st2e} object, which contains a matrix representing variable selection
#' ensemble and a numerical vector of objective scores corresponding to each row of the ensemble that represents one \code{st2} selection.
#' Importance scores are calculated as selection frequencies weighted by objective scores.
#'
#' @param st2e_obj A \code{st2e} object returned by \code{\link[StableMate]{st2e}}.
#' @param st2e_obj2 A \code{st2e} object whose ensemble is the predictor pool of \code{st2e_obj}. If provided, the
#' importance score calculation will be adjusted accordingly.
#' @param prune Logical; If TRUE, weight selections by ranking quantiles of their objective scores when calculating importance scores.
#' @param prefilt_scores A numerical vector that contains some weighting measurements of predictor pools. Predictor pools will be
#' weighted by ranking quantiles of \code{prefilt_scores} when calculating importance scores.
#' When \code{st2e_obj2} is given, \code{prefilt_scores} will be set as the objective scores of \code{st2e_obj2}. By default is set to
#' a vector of ones.
#' @param scale_psd_imp Logical; If TRUE, scale the importance score of the pseudo-predictor by importance scores calculated on
#' predictor pools. Used only when the predictor pools of \code{st2e_obj} are pre-filtered by methods other than  \code{\link[StableMate]{st2e}}.
#' @param per_pred_evl Logical; If TRUE, each predictor will be benchmarked against the pseudo-predictor selections whose
#' predictor pools contains the query predictor. Otherwise, all predictors are benchmarked against the same set of pseudo-predictor selections,
#' which are usually all the selections in the ensemble if the pseudo-predictor is included in every predictor pool.
#' @param B Numeircal; The number of bootstrapping iteration for calculating importance scores.
#'
#' @return A \code{st2e} object with updated importance scores (stored in the \code{imp_scores} of the \code{st2e} object).
#'
#' @export
#'
#' @seealso \code{\link[StableMate]{st2e}()}
#' @name calc_imp
#' @rdname calc_imp
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
