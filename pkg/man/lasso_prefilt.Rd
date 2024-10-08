% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lasso_prefilt.R
\name{lasso_prefilt}
\alias{lasso_prefilt}
\title{Build pre-filtering ensemble by a random Lasso procedure.}
\usage{
lasso_prefilt(
  Y,
  X,
  p,
  env = NULL,
  K = 100,
  ncore = NULL,
  drop_pred = TRUE,
  par_method = c("SNOW", "MC", "MPI"),
  verbose = TRUE,
  chunk_size = ceiling(K/ncore),
  group_by = c("env", "Yenv", "none"),
  ...
)
}
\arguments{
\item{Y}{A response vector or matrix depending on the objective function. It should be a vector if the function is used with
StableMate's default objective.}

\item{X}{A predictor matrix with rows representing samples and columns representing predictors. The columns must be named.}

\item{p}{Numerical; Pre-filter size for each random Lasso procedure.}

\item{K}{Numerical; Ensemble size.}

\item{ncore}{Numerical; Numerical; If greater than 0. Parallel computing is enabled.}

\item{drop_pred}{Logical; If true, remove predictors from X that are not in the predictor pool to save space and time.}

\item{verbose}{Logical; If true, multi-core verbose will be printed.}

\item{chunk_size}{Numerical; The size of task chunks (How many st2 runs are performed one each task).}

\item{group_by}{Either 'env', 'Yenv' or 'none' to weight samples by the inverse sizes of their environments,
their Y class-environments combined if binomial regression, or no weighting respectively.}

\item{...}{Arguments passed to \code{\link[glmnet]{glmnet}}.}

\item{par_method;}{Parallel computing method. SNOW is preferred on Windows local machines, MC is preferred on non-Windows local machines.
MPI is preferred on distributed computing cloud.}
}
\description{
Build pre-filtering ensemble by a random Lasso procedure that is similar to
\href{https://doi.org/10.1111/j.1467-9868.2010.00740.x}{Stability selection, Meinshausen et al. (2010, The Journal of the Royal Statistical Society, Series B)})
For each st2 run, we randomly sample one half of samples to select the top \code{p} predictors with Lasso.
The pre-filtered predictors is then set as the predictor pool. The advantages are of two folds.
First, across the different re-sampling runs, the top \code{p} predictors are expected to differ thus
enabling to cover a large and diverse range of predictors in our overall search.
Second, we improve the stability of Lasso pre-screening when subjected to sample perturbation.
}
