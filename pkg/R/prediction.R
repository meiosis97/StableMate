# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Built-in \code{st2e} prediction functions.
#'
#' Given a regression model that is returned by \code{\link[StableMate]{reg_fun}} provided to \code{\link[StableMate]{st2e}},
#' and given new predictor data, makes predictions on the response.
#'
#' @param model A regression model that is returned by \code{\link[StableMate]{reg_fun}} provided to \code{\link[StableMate]{st2e}}.
#' @param X New predictor data that can be recognized by \code{model}
#'
#' @return Predicted Y.
#'
#' @export
#'
#' @seealso \code{\link[StableMate]{reg_fun}}
#' @name pred_fun
#' @rdname pred_fun
pred_ols <- function(model, X){
  model[1] + X %*% model[colnames(X)]
}


#' Logistic regression
#'
#' @export
#' @rdname pred_fun
pred_logit <- function(model, X){
  ilogit <- function(x) exp(x)/(1+exp(x))
  ilogit(model[1] + X %*% model[colnames(X)])
}
