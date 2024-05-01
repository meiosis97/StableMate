# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' StableMate regression and stability analysis
#'
#' Perform stabilized regression (\href{https://doi.org/10.1214/21-AOAS1487}{Pfister et al. (2019, The Annals of Applied Statistics)})
#' by stochastic stepwise variable selection. Many controls in this function are shared by the \code{\link[StableMate]{st2e}}.
#'
#' @param Y A response vector or matrix depending on the objective function. It should be a vector if the function is used with
#' StableMate's default objective.
#' @param X A predictor matrix with rows representing samples and columns representing predictors. The columns must be named.
#' @param env A character vector indicates sample environments. Should be of the same length as the number of rows of \code{X}.
#' @param K Numerical; Ensemble size.
#' @param max_size Numerical; Any step of selection will not consider a predictor set that is larger than \code{max_size}.
#' @param alpha Numerical; A Parameter of beta distribution for sampling proposal size.
#' @param beta Numerical; A Parameter of beta distribution for sampling proposal size.
#' @param lambda Numerical; A scaling factor that controls the size of binomial distribution for sampling proposal size.. Should be a value between 0 and 1.
#' If 1, then the proposal size can be as large as the size of the maximum possible proposal.
#' @param a Numerical; A scaling factor that controls the number of proposal to sample, which contains a term where \code{a} is multiplied by the size of the maximum possible proposal
#' to the power of \code{t}.
#' @param b Numerical; A scaling factor that controls the number of proposal to sample, which contains a term where \code{b} is multiplied by a quality that measures the number of
#' all possible unique proposal that can be made.
#' @param t Numerical; A scaling factor that controls the number of proposal to sample, which contains a term where the size of the maximum possible proposal is raised to the power of
#' \code{t} and multiplied by \code{a}.
#' @param pred_st2_ctrl Controls for \code{\link[StableMate]{st2e}} that selects most predictive features of the response. See \code{\link[StableMate]{st2e}} for the detailed description of controls.
#' \code{K}, \code{ret_imp}, \code{calc_imp_ctrl}, \code{drop_pred} will be ignored if provided.
#' @param stab_st2_ctrl Controls for \code{\link[StableMate]{st2e}} that selects most stable preidctors of the response. See \code{\link[StableMate]{st2e}} for the detailed description of controls.
#' \code{K}, \code{ret_imp}, \code{calc_imp_ctrl}, \code{drop_pred} will be ignored if provided. By default, stable predictors will be searched within predictive features.
#' @param pred_calc_imp_ctrl Controls for \code{\link[StableMate]{calc_imp}} applied on the predictive ensemble for calculating predictivity scores.
#' See \code{\link[StableMate]{calc_imp}} for the detailed description of controls.
#' @param stab_calc_imp_ctrl Controls for \code{\link[StableMate]{calc_imp}} applied on the stable and predictive ensemble for calculating stability scores.
#' See \code{\link[StableMate]{calc_imp}} for the detailed description of controls.  By defaults, selections in the stable and predictive ensemble
#' are considered as the subsets of selections in the predictive ensemble, and stability scores are adjusted accordingly.
#' @param do_switch Logical; If true, switch steps that swap predictors that are in and out the model are allowed.
#' @param ret_mod Logical; If true, return the prediction ensemble
#' @param ret_imp Logical; If true, return the importance score.
#' @param ncore Numerical; Numerical; If greater than 0. Parallel computing is enabled.
#' @param par_method; Parallel computing method. SNOW is preferred on Windows local machines, MC is preferred on non-Windows local machines.
#' MPI is preferred on distributed computing cloud.
#' @param verbose Logical; If true, multi-core verbose will be printed.
#' @param chunk_size Numerical; The size of task chunks (How many st2 runs are performed one each task).
#' @param fun_export A vector of names of functions to be imported by \code{\link[foreach]{foreach}} environment.
#' @param pkg_export A vector of names of packages to be imported by\code{\link[foreach]{foreach}} environment.
#' @param drop_pred Logical; If true, remove predictors from X that are not in the predictor pool to save space and time.
#'
#' @return A \code{stablemate} object.
#'
#' @export
#' @rdname stablemate
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
      if(!is.null(attributes(pred_st2_ctrl$pred_pool)$stbm_lasso) & is.null(pred_calc_imp_ctrl$scale)){
        pred_calc_imp_ctrl$scale <- TRUE

      }

      # A list containing all the default controls for calculating importance scores
      defaults <- list(prune = T, prefilt_scores = NULL, scale = F, pooled = T, B = 5000)

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
      if(class(pred_st2_ctrl$pred_pool) == "lasso_ensemble" & is.null(pred_calc_imp_ctrl$scale)){
        pred_calc_imp_ctrl$scale <- TRUE

      }

      # A list containing all the default controls for calculating importance scores measuring predictivity
      defaults <- list(prune = T, prefilt_scores = NULL, scale = F, pooled = T, B = 5000)

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
      if(class(stab_st2_ctrl$pred_pool) == "lasso_ensemble" & is.null(stab_calc_imp_ctrl$scale)){
        stab_calc_imp_ctrl$scale <- TRUE

      }

      # A list containing all the default controls for calculating the importance scores measuring stability of predictors
      defaults <- list(prune = T, prefilt_scores = NULL, scale = F, pooled = T, B = 5000)

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
