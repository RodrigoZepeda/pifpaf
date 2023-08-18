#' @title Potential Impact Fraction via bootstrap for survey data
#'
#' @description Estimates the **potential impact fraction**, `pif`, for
#' **individual-level exposure** data (and covariates), `X`, from a
#' cross-sectional survey. Exposure is assumed to be associated with
#' a **relative risk** function, `rr`, with parameter `theta`. A **counterfactual**
#' scenario as a function of the exposure `cft` is assumed.
#'
#' The **potential impact fraction** is defined \insertCite{chan2023nonparametric;textual}{pifpaf}:
#' \loadmathjax
#' \mjdeqn{
#'  \text{PIF} = \dfrac{\mathbb{E}\Big[RR(X;\theta)\Big] -
#'    \mathbb{E}\Big[RR\big(\text{cft}(X);\theta\big)\Big]}{\mathbb{E}\Big[RR(X;\theta)\Big]}
#' }{pif = (mean(rr) - mean(cft))/mean(rr)}
#'
#' where:
#' * \mjseqn{X} denotes the individual-level matrix of exposure and covariates,
#' * \mjseqn{\theta} represents additional parameters of the relative risk function,
#' * \mjseqn{RR(X,\theta)} denotes the relative risk of exposure (and covariates) at level
#' \mjseqn{X} given parameters \mjseqn{\theta},
#' * \mjseqn{cft(X)} denotes the counterfactual function applied to exposure and covariates,
#' * \mjseqn{{\mathbb{E}\Big[RR(X;\theta)\Big]}} and
#' \mjseqn{{\mathbb{E}\Big[RR(\text{cft}(X);\theta)\Big]}} denote the population average relative
#' risk under current (observed) conditions and the relative risk under the counterfactual scenario.
#' * \mjseqn{\text{PIF}} represents the potential impact fraction.
#'
#' @param design (`survey.design`/`svyrep.design`) survey design structure from
#' the survey package obtained from [survey::svydesign()] or a survey replicates design
#' obtained either from [survey::svrepdesign()] or [svrep::as_bootstrap_design()]. Contains
#' individual level exposure and covariates.
#'
#' @param theta (`vector`/`double`) parameters of the relative risk function `rr`.
#'
#' @param rr (`function`) relative risk function with two parameters: `X` and `theta` (in that order).
#'
#' @param cft (`function`) counterfactual function with one parameter: `X`.
#'
#' @param n_samples (`double`) number of bootstrap samples
#'
#' @param additional_theta_arguments any additional information on `theta` utilized
#' for obtaining bootstrap samples from the paramter. Options are:
#'  + (`double`) the **variance** of `theta` if `theta` is
#'  one dimensional and asymptotical normality is assumed (default).
#'  + (`vector`) the **variances** of each entry of `theta` if `theta` is
#'  n-dimensional and its entries are uncorrelated and asymptotical normality is assumed (default).
#'  + (`matrix`) the **variance-covariance** matrix of `theta` if `theta` is
#'  n-dimensional and its entries are correlated and asymptotical normality is assumed (default).
#'  + any argument to pass to `theta_distribution` to simulate samples from `theta`
#'  if `theta` is not assumed to be asymptotically normally distributed.
#'
#' **Optional**
#'
#' @param theta_distribution (`function`) random number generator that follows
#' the distribution of the estimator `theta`. By default, `theta` is assumed to be asymptotically
#' normal and thus `theta_distribution` is set to [mvtnorm::rmvnorm()] with
#' variance given by `additional_theta_arguments`.
#'
#' @param uncertainty_interval_type (`string`) either `'wald'` (recommended) or `'percentile'` for
#' uncertainty intervals based upon the bootstrap's percentile or the Wald-type approximation
#' (see confidence interval section)
#'
#' @param parallel (`boolean`) Flag indicating whether or not to do computations in parallel.
#' Default is `TRUE`.
#'
#' @param num_cores (`int`) Number of cores; defaults to [parallel::detectCores()].
#'
#' @param confidence_level (`double`) Confidence level for the uncertainty interval.
#' Defaults to `0.95`
#'
#' @section Confidence Intervals:
#'
#' Confidence intervals are estimated via the two methods specified
#' by \insertCite{arnab2017survey;textual}{pifpaf}.
#'
#' 1. **Wald-type confidence intervals** (which are more precise) are of the form:
#' \mjdeqn{
#' \widehat{\text{PIF}} \pm t_{1 - \alpha/2}
#' \sqrt{\widehat{\text{Var}}_{\text{B}}\big[\widehat{\text{PIF}}\big]}
#' }{pif +- qt(1 - alpha/2)*sqrt(var(pif))}
#' where \mjseqn{t_{1 - \alpha/2}} stands for the percent points at level \mjseqn{1 - \alpha/2}
#' of Student's t-distribution, and
#' \mjseqn{\widehat{\text{Var}}_{\text{B}}\big[\widehat{\text{PIF}}\big]}
#' represents the estimated variance (via bootstrap) of the potential impact fraction estimator.
#'
#' 2. **Percentile confidence intervals** (less precise) are of the form:
#'\mjdeqn{
#' \Big[\widehat{\text{PIF}}_{\text{B},\alpha/2}, \widehat{\text{PIF}}_{\text{B},1-\alpha/2}\Big]
#'}{quantile(pif, c(alpha/2, 1-alpha/2))}
#' where \mjseqn{\widehat{\text{PIF}}_{\text{B},k}} represents the kth sample quantile
#' from the bootstrap simulation of the potential impact fraction estimator.
#'
#' @references \insertAllCited{}
#'
#' @seealso [paf()], [pif()], [pif_approximate()], [pif_delta_method()].
#'
#' @keywords internal, pif, impact-fraction
#' @export
#' @md

pif_survey_bootstrap <- function(){

  #Get the weights
  #https://stackoverflow.com/questions/73627746/how-do-i-get-the-weights-from-a-survey-design-object-in-r
  weights(rep_survey_design, type = "analysis")
}
