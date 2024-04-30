# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Built-in \code{st2e} objective functions
#'
#' Evaluate the objective score of a given set of predictors (proposals). \code{\link[StableMate]{st2e}} algorithm
#' searches for the best of predictors defined by \code{\link[StableMate]{sel_fun}}. The customized objective function should
#' have same input as described below, and output a objective score
#'
#'
#' @param Y A response vector or matrix depending on the objective function. It should be a vector if the function is used with
#' StableMate's default objective.
#' @param X A predictor matrix with rows representing samples and columns representing predictors. The columns must be named.
#' @param env A character vector indicates sample environments. Should be of the same length as the number of rows of \code{X}.
#' @param pred_in A indicator vector of the same length as the number of columns of \code{X}. Ones in the vector indicates the
#' corresponding predictors matched in column positions to be evaluated.
#'
#' @return A objective score of a minimization problem.
#'
#' @export
#'
#' @seealso \code{\link[StableMate]{st2e}()}
#' @name obj_fun
#' @rdname obj_fun
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


#' OLS leave-one-environment-out cross validation sum of squares of errors.
#'
#' @export
#' @rdname obj_fun
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


#' Logistic regression Bayesian information criteria.
#'
#' @export
#' @rdname obj_fun
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


#' Logistic regression leave-one-environment-out cross validation sum of deviance.
#'
#' @export
#' @rdname obj_fun
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
