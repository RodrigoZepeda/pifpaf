% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pifpaf_class.R
\name{coef.pif_class}
\alias{coef.pif_class}
\title{Extract point estimate of potential impact fraction or population attributable fraction}
\usage{
\method{coef}{pif_class}(x, ...)
}
\arguments{
\item{x}{A \code{pif_class} object}

\item{...}{additional arguments for \code{\link[base:mean]{base::mean()}}}
}
\description{
Method for obtaining the point estimate of a \code{pif_class} object (i.e. the
point-estimate of the potential impact fraction or population attributable fraction via
the bootstrap).
}
\examples{
#Example 1
data(ensanut)
options(survey.lonely.psu = "adjust")
design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
rr <- function(X, theta) {
  exp(-2 +
    theta[1] * X[, "age"] + theta[2] * X[, "systolic_blood_pressure"] / 100)
}
cft <- function(X) {
  X[, "systolic_blood_pressure"] <- X[, "systolic_blood_pressure"] - 5
  return(X)
}
pifsim <- pif(design,
  theta = log(c(1.05, 1.38)), rr, cft,
  additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10,
)
coef(pifsim)
}
