% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/st2e_util.R
\name{print.st2e}
\alias{print.st2e}
\title{print function for \code{st2e} object.}
\usage{
\method{print}{st2e}(st2e_obj, imp_type = c("conditional", "joint"), sigthresh = 1 - exp(-5))
}
\arguments{
\item{st2e_obj}{A \code{st2e} object returned by \code{\link[StableMate]{st2e}}.}

\item{imp_type}{Either 'conditional' or 'joint' to make selections based on conditional importance scores
or joint importance scores.}

\item{sigthresh}{Numerical; Significance threshold.}
}
\description{
print function for \code{st2e} object.
}
