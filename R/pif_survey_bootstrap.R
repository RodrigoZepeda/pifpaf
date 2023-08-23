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
#' @inheritParams pif
#' @inheritSection pif Confidence Intervals
#' @inheritSection pif Additional parallelization options
#'
#' @references \insertAllCited{}
#'
#' @keywords internal
#' @importFrom parallel detectCores
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @md

pif_survey_bootstrap <- function(design,
                                 theta,
                                 rr,
                                 cft = NULL,
                                 additional_theta_arguments,
                                 n_bootstrap_samples = NULL,
                                 theta_distribution = "default",
                                 uncertainty_interval_type = c("wald", "percentile"),
                                 parallel = TRUE,
                                 num_cores = 1,
                                 confidence_level = 0.95,
                                 return_replicates = FALSE,
                                 is_paf = FALSE,
                                 ...) {

  # Validators for almost all of the parameters
  # theta can be literally anything so we don't validate
  rr                  <- validate_rr(rr = rr)
  is_paf              <- validate_is_paf(is_paf = is_paf)
  cft                 <- validate_cft(cft = cft, is_paf = is_paf)
  confidence_level    <- validate_confidence_level(confidence_level = confidence_level)
  num_cores           <- validate_number_of_cores(num_cores = num_cores)
  `%dofun%`           <- validate_parallel_setup(parallel = parallel, num_cores = num_cores)
  return_replicates   <- validate_return_replicates(return_replicates = return_replicates)
  theta_args          <- validate_theta_arguments(theta_distribution = theta_distribution,
                            additional_theta_arguments = additional_theta_arguments, theta = theta)
  design             <- validate_survey_design(design = design,
                            n_bootstrap_samples = n_bootstrap_samples, ...)
  theta_distribution <- validate_theta_distribution(theta_distribution = theta_distribution,
                            additional_theta_arguments = theta_args)


  # Make cluster
  if (!is.null(num_cores) & num_cores > 1) {
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
  }

  # Get the weights and data from the survey
  # https://stackoverflow.com/questions/73627746/how-do-i-get-the-weights-from-a-survey-design-object-in-r
  data <- design[["variables"]]
  weight_matrix <- stats::weights(design, type = "analysis")

  # Get the simulated thetas
  theta_args <- append(theta_args, list("n" = ncol(weight_matrix)))
  sim_theta <- do.call(theta_distribution, args = theta_args)

  # Loop through the foreach
  i <- 0 # Declare the global
  pif_replicates <- foreach::foreach(
    i = 1:ncol(weight_matrix),
    .combine = c,
    .inorder = FALSE,
    .multicombine = TRUE,
    .noexport = c("design")
  ) %dofun% {
    # Obtain the pif
    internal_pif_data_frame(
      df = data,
      theta = sim_theta[i, ],
      rr = rr,
      cft = cft,
      weights = weight_matrix[, i],
      is_paf = is_paf
    )
  }

  if (!is.null(num_cores) & num_cores > 1) {
    doParallel::stopImplicitCluster(cl)
  }

  # Get point estimate
  point_estimate <- mean(pif_replicates)
  std_pif <- sd(pif_replicates)

  # Get alpha of Uncertainty interval
  alpha_ui <- 1 - confidence_level

  if (uncertainty_interval_type[1] == "percentile") {
    interval <- quantile(pif_replicates, c(alpha_ui / 2, 1 - alpha_ui / 2))
    names(interval) <- c("Lower", "Upper")
  } else if (uncertainty_interval_type[1] == "wald") {
    interval <- c(
      "Lower" = point_estimate - qt(1 - alpha_ui / 2, df = ncol(weight_matrix)) * std_pif,
      "Upper" = point_estimate + qt(1 - alpha_ui / 2, df = ncol(weight_matrix)) * std_pif
    )
  }

  pif_simulations <- list(
    "Point" = point_estimate,
    "Interval" = interval,
    "Confidence" = confidence_level,
    "Standard_Deviation" = std_pif,
    "Uncertainty interval type" = uncertainty_interval_type[1]
  )

  if (return_replicates) {
    pif_simulations <- append(pif_simulations, list("Replicates" = pif_replicates))
  }

  return(pif_simulations)
}


#' @title Pif from dataframe
#' @keywords internal
#' @note In previous version of the package this function was called `pif.empirical`.
internal_pif_data_frame <- function(df, theta, rr, cft,
                                    weights = NULL,
                                    is_paf = FALSE) {
  # Estimate weighted sums
  if (is.null(weights)) {
    mux <- mean(as.matrix(rr(df, theta)))
  } else {
    mux <- weighted.mean(as.matrix(rr(df, theta)), weights)
  }


  if (is_paf) {
    mucft <- 1
  } else {
    if (is.null(weights)) {
      mucft <- mean(as.matrix(rr(cft(df), theta)))
    } else {
      mucft <- weighted.mean(as.matrix(rr(cft(df), theta)), weights)
    }
  }

  # Calculate PIF
  if (is.infinite(mux)) {
    cli::cli_alert_warning("Expected value of Relative Risk is not finite")
  }

  if (is.infinite(mucft)) {
    cli::cli_alert_warning("Expected value of Relative Risk under counterfactual is not finite")
  }

  pif <- 1 - mucft / mux

  return(pif)
}

#' @title Potential Impact Fraction
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
#' @param design (`survey.design`, `data.frame`,`tibble`, or `svyrep.design`) survey data
#' structure. If data comes from a survey set the `design` with [survey::svydesign()]. It can
#' also support a [survey::svrepdesign()] design if your survey comes with replicates.
#' Finally, the model can also accommodate a `data.frame` or `tibble` with weights assuming
#' simple random sampling without replacement.
#'
#' @param theta (`vector`/`double`) parameters of the relative risk function `rr`.
#'
#' @param rr (`function`) relative risk function with two parameters: a `data.frame` called
#' `X` containing the individual-level exposure and covariates, and `theta` (in that order).
#'
#' @param cft (`function`) counterfactual function that takes a `data.frame`, `X`, of
#' individual-level exposure and covariates and returns a new `data.frame` of individual-level
#' counterfactual exposure and covariates.
#'
#' @param additional_theta_arguments any additional information on `theta` utilized
#' for obtaining bootstrap samples from the paramter. Options are:
#'  + (`double`) the **variance** of `theta` if `theta` is
#'  one dimensional and asymptotical normality is assumed (default).
#'  + (`vector`) the **variances** of each entry of `theta` if `theta` is
#'  n-dimensional and its entries are uncorrelated and asymptotical normality is assumed (default).
#'  + (`matrix`) the **variance-covariance** matrix of `theta` if `theta` is
#'  n-dimensional and its entries are correlated and asymptotical normality is assumed (default).
#'  + any list of arguments to pass via [base::do.call()] to `theta_distribution` to
#'  simulate samples from `theta` if `theta` is not assumed to be asymptotically normally
#'  distributed.
#'
#' **Optional**
#'
#' @param n_bootstrap_samples (`double`) number of bootstrap samples. If a `svyrep.design` is passed
#' as an argument, then `n_bootstrap_samples` represents the number of number of replicates in
#' the design.
#'
#' @param theta_distribution (`function`) random number generator that follows
#' the distribution of the estimator `theta`. By default, `theta` is assumed to be asymptotically
#' normal and thus `theta_distribution` is set to [mvtnorm::rmvnorm()] with
#' variance given by `additional_theta_arguments`. The number of simulations for the
#' `theta_distribution` function must be parametrized by a parameter of name `n`.
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
#' Defaults to `0.95`.
#'
#' @param is_paf (`boolean`) Whether the function being estimated is the Population
#' Attributable Fraction (`is_paf = TRUE`) or the Potential Impact Fraction
#' (`is_paf = FALSE`)
#'
#' @param return_replicates (`boolean`) Whether to return the simulated impact fractions
#' or not.
#'
#' @param weights  (`vector`) If you are not following the recommended version and use a
#' `svydesign` object for the design you can still use `weights` to associate
#' weights to your estimation. Beware that it might not give accurate estimations of the variance
#' nor the uncertainty intervals.
#'
#' @param ... Additional parameters for [svrep::as_bootstrap_design()].
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
#' \mjdeqn{
#' \Big[\widehat{\text{PIF}}_{\text{B},\alpha/2}, \widehat{\text{PIF}}_{\text{B},1-\alpha/2}\Big]
#' }{quantile(pif, c(alpha/2, 1-alpha/2))}
#' where \mjseqn{\widehat{\text{PIF}}_{\text{B},k}} represents the kth sample quantile
#' from the bootstrap simulation of the potential impact fraction estimator.
#'
#' @section Additional parallelization options:
#' By default the function uses [foreach::foreach] to parallelize creating a cluster
#' via [doParallel::registerDoParallel]. For finner parallelization control you can do:
#' ```
#' #Create the cluster outside using whatever you want
#' myCluster <- parallel::makeCluster(3, type = "PSOCK")
#'
#' #Register the cluster
#' doParallel::registerDoParallel(myCluster)
#'
#' #Call the function using parallel = TRUE to use %dopar% and num_cores = NULL to avoid setting
#' pif(..., parallel = TRUE, num_cores = NULL)
#'
#' #Stop the cluster
#' doParallel::stopCluster(myCluster)
#' ```
#' @examples
#' # Use the ensanut dataset
#' data(ensanut)
#'
#' # EXAMPLE 1
#' # Setup the survey design
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#' rr <- function(X, theta) {
#'   exp(-2 +
#'     theta[1] * X[, "age"] + theta[2] * X[, "systolic_blood_pressure"] / 100)
#' }
#' cft <- function(X) {
#'   X[, "systolic_blood_pressure"] <- X[, "systolic_blood_pressure"] - 5
#'   return(X)
#' }
#' pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#'
#' # EXAMPLE 2
#' # Now do the same but using a replicate design
#' options(survey.lonely.psu = "adjust")
#' rep_design <- svrep::as_bootstrap_design(design, replicates = 10)
#' pif(rep_design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03)
#' )
#' @references \insertAllCited{}
#' @seealso [paf()]
#' @keywords pif impact-fraction
#' @export
#' @importFrom parallel detectCores
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @md

pif <- function(design,
                theta,
                rr,
                cft = NULL,
                additional_theta_arguments,
                n_bootstrap_samples = NULL,
                theta_distribution = "default",
                uncertainty_interval_type = c("wald", "percentile"),
                parallel = TRUE,
                num_cores = 1,
                confidence_level = 0.95,
                return_replicates = FALSE,
                is_paf = FALSE,
                weights = NULL,
                ...) {


  #If is data frame or tibble create a survey design that does simple random sampling
  if (inherits(design, "data.frame")){
    design <- survey::svydesign(id =~1, data = design, weights = weights)
  }

  return(
    pif_survey_bootstrap(design = design,
                         theta = theta,
                         rr = rr,
                         cft = cft,
                         additional_theta_arguments = additional_theta_arguments,
                         n_bootstrap_samples = n_bootstrap_samples,
                         theta_distribution = theta_distribution,
                         uncertainty_interval_type = uncertainty_interval_type,
                         parallel = parallel, num_cores = num_cores,
                         confidence_level = confidence_level,
                         return_replicates = return_replicates,
                         is_paf = is_paf, ...)
    )
}

#' @title Population Attributable Fraction
#'
#' @description Estimates the **population attributable fraction**, `paf`, for
#' **individual-level exposure** data (and covariates), `X`, from a
#' cross-sectional survey. Exposure is assumed to be associated with
#' a **relative risk** function, `rr`, with parameter `theta`. A **counterfactual**
#' scenario as a function of the exposure `cft` is assumed.
#'
#' The **population attributable fraction** is defined
#' \insertCite{chan2023nonparametric;textual}{pifpaf} as:
#' \loadmathjax
#' \mjdeqn{
#'  \text{PAF} = \dfrac{\mathbb{E}\Big[RR(X;\theta)\Big] - 1}{\mathbb{E}\Big[RR(X;\theta)\Big]}
#' }{pif = (mean(rr) - 1)/mean(rr)}
#'
#' where:
#' * \mjseqn{X} denotes the individual-level matrix of exposure and covariates,
#' * \mjseqn{\theta} represents additional parameters of the relative risk function,
#' * \mjseqn{RR(X,\theta)} denotes the relative risk of exposure (and covariates) at level
#' \mjseqn{X} given parameters \mjseqn{\theta},
#' * \mjseqn{{\mathbb{E}\Big[RR(X;\theta)\Big]}} represents the average relative risk
#' in the population
#' * \mjseqn{1} represents the relative risk under the theoretical minimum risk scenario,
#' * \mjseqn{\text{PAF}} represents the population attributable fraction
#'
#' @inheritParams pif
#'
#' @section Confidence Intervals:
#'
#' Confidence intervals are estimated via the two methods specified
#' by \insertCite{arnab2017survey;textual}{pifpaf}.
#'
#' 1. **Wald-type confidence intervals** (which are more precise) are of the form:
#' \mjdeqn{
#' \widehat{\text{PAF}} \pm t_{1 - \alpha/2}
#' \sqrt{\widehat{\text{Var}}_{\text{B}}\big[\widehat{\text{PIF}}\big]}
#' }{pif +- qt(1 - alpha/2)*sqrt(var(PAF))}
#' where \mjseqn{t_{1 - \alpha/2}} stands for the percent points at level \mjseqn{1 - \alpha/2}
#' of Student's t-distribution, and
#' \mjseqn{\widehat{\text{Var}}_{\text{B}}\big[\widehat{\text{PIF}}\big]}
#' represents the estimated variance (via bootstrap) of the potential impact fraction estimator.
#'
#' 2. **Percentile confidence intervals** (less precise) are of the form:
#' \mjdeqn{
#' \Big[\widehat{\text{PAF}}_{\text{B},\alpha/2}, \widehat{\text{PAF}}_{\text{B},1-\alpha/2}\Big]
#' }{quantile(pif, c(alpha/2, 1-alpha/2))}
#' where \mjseqn{\widehat{\text{PAF}}_{\text{B},k}} represents the kth sample quantile
#' from the bootstrap simulation of the potential impact fraction estimator.
#'
#' @section Additional parallelization options:
#' By default the function uses [foreach::foreach] to parallelize creating a cluster
#' via [doParallel::registerDoParallel]. For finner parallelization control you can do:
#' ```
#' #Create the cluster outside using whatever you want
#' myCluster <- parallel::makeCluster(3, type = "PSOCK")
#'
#' #Register the cluster
#' doParallel::registerDoParallel(myCluster)
#'
#' #Call the function using parallel = TRUE to use %dopar% and num_cores = NULL to avoid setting
#' paf(..., parallel = TRUE, num_cores = NULL)
#'
#' #Stop the cluster
#' doParallel::stopCluster(myCluster)
#' ```
#' @examples
#' # Use the ensanut dataset
#' data(ensanut)
#'
#' # EXAMPLE 1
#' # Setup the survey design
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#' rr <- function(X, theta) {
#'   exp(-2 +
#'     theta[1] * X[, "age"] + theta[2] * X[, "systolic_blood_pressure"] / 100)
#' }
#' paf(design,
#'   theta = log(c(1.05, 1.38)), rr,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#'
#' # EXAMPLE 2
#' # Now do the same but using a replicate design
#' options(survey.lonely.psu = "adjust")
#' rep_design <- svrep::as_bootstrap_design(design, replicates = 10)
#' paf(rep_design,
#'   theta = log(c(1.05, 1.38)), rr,
#'   additional_theta_arguments = c(0.01, 0.03)
#' )
#' @references \insertAllCited{}
#' @seealso [pif()]
#' @keywords paf attributable-fraction
#' @export
#' @importFrom parallel detectCores
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @md

paf <- function(design,
                theta,
                rr,
                additional_theta_arguments,
                n_bootstrap_samples = NULL,
                theta_distribution = "default",
                uncertainty_interval_type = c("wald", "percentile"),
                parallel = TRUE,
                num_cores = 1,
                confidence_level = 0.95,
                return_replicates = FALSE,
                weights = NULL,
                ...) {


  return(
    pif(design = design,
        theta = theta,
        rr = rr,
        additional_theta_arguments = additional_theta_arguments,
        n_bootstrap_samples = n_bootstrap_samples,
        theta_distribution = theta_distribution,
        uncertainty_interval_type = uncertainty_interval_type,
        parallel = parallel, num_cores = num_cores,
        confidence_level = confidence_level,
        return_replicates = return_replicates,
        is_paf = TRUE,
        weights = weights,
         ...)
  )
}


