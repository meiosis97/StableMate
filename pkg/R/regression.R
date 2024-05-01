# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Built-in \code{st2e} regression functions.
#'
#' Build regression models based on the optimal set of predictors (optimal proposal) after \code{\link[StableMate]{st2e}} search.
#' The customized regression function should have same input as described below, and output a regression model that can be
#' used to make predictions.
#'
#' @param Y A response vector or matrix depending on the objective function. It should be a vector if the function is used with
#' StableMate's default objective.
#' @param X A predictor matrix with rows representing samples and columns representing predictors. The columns must be named.
#' @param env A character vector indicates sample environments. Should be of the same length as the number of rows of \code{X}.
#' @param pred_in A indicator vector of the same length as the number of columns of \code{X}. Ones in the vector indicates the
#' corresponding predictors matched in column positions to be fitted.
#'
#' @return A regression model.
#'
#' @export
#'
#' @seealso \code{\link[StableMate]{st2e}()}
#' @name reg_fun
#' @rdname reg_fun
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


#' Logistic regression
#'
#' @export
#' @rdname reg_fun
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
