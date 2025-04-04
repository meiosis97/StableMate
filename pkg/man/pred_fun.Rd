% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.R
\name{pred_fun}
\alias{pred_fun}
\alias{pred_ols}
\alias{pred_logit}
\title{Built-in \code{st2e} prediction functions.}
\usage{
pred_ols(model, X)

pred_logit(model, X)
}
\arguments{
\item{model}{A regression model that is returned by \code{\link[StableMate]{reg_fun}} provided to \code{\link[StableMate]{st2e}}.}

\item{X}{New predictor data that can be recognized by \code{model}}
}
\value{
Predicted Y.
}
\description{
Given a regression model that is returned by \code{\link[StableMate]{reg_fun}} provided to \code{\link[StableMate]{st2e}},
and given new predictor data, makes predictions on the response.
}
\seealso{
\code{\link[StableMate]{reg_fun}}
}
