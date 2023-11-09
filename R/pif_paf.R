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
#' @inheritSection pif Additional parallelization options
#'
#' @references \insertAllCited{}
#'
#' @keywords internal
#' @importFrom doFuture  %dofuture%
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @seealso [pif()] [paf()]
#' @return A `pif_class` object estimating the potential impact fraction for individual-level data
#' from `design` using relative risk function `rr` and counterfactual `cft`.
#' @md

pif_survey_bootstrap <- function(design,
                                 theta,
                                 rr,
                                 cft = NULL,
                                 additional_theta_arguments,
                                 n_bootstrap_samples = NULL,
                                 theta_distribution = "default",
                                 is_paf = FALSE,
                                 .options.future = list(seed = TRUE),
                                 ...) {

  # Validators for almost all of the parameters
  # theta can be literally anything so we don't validate
  rr                  <- validate_rr(rr = rr)
  is_paf              <- validate_is_paf(is_paf = is_paf)
  cft                 <- validate_cft(cft = cft, is_paf = is_paf)
  theta_args          <- validate_theta_arguments(theta_distribution = theta_distribution,
                            additional_theta_arguments = additional_theta_arguments, theta = theta)
  design             <- validate_survey_design(design = design,
                            n_bootstrap_samples = n_bootstrap_samples, ...)
  theta_distribution <- validate_theta_distribution(theta_distribution = theta_distribution,
                            additional_theta_arguments = theta_args)


  # Get the weights and data from the survey
  # https://stackoverflow.com/questions/73627746/how-do-i-get-the-weights-from-a-survey-design-object-in-r
  data <- design[["variables"]]
  weight_matrix <- stats::weights(design, type = "analysis")

  # Get the simulated thetas
  theta_args <- append(theta_args, list("n" = ncol(weight_matrix)))
  sim_theta  <- do.call(theta_distribution, args = theta_args)

  # Loop through the foreach
  boot_i <- 0; rr_i <- 0; cft_i <- 0  # Declare the global for CRAN

  cli::cli_progress_bar("Simulating", total = length(rr)*length(cft))
  #Loop through each relative risk
  pif_values <- foreach::foreach(rr_i = 1:length(rr), .combine = rbind) %do% {
    #Loop through each counterfactual
    foreach::foreach(cft_i = 1:length(cft), .combine = rbind,
                     .options.future = .options.future) %dofuture% {
      #Loop through each replicate
      pif_replicates <- foreach::foreach(boot_i = 1:ncol(weight_matrix), .combine = rbind,
                                         .options.future = .options.future,
                                         .inorder = FALSE, .multicombine = TRUE) %dofuture% {
        #Obtain the potential impact fraction
        internal_pif_data_frame(df = data, theta = sim_theta[boot_i, ], rr = rr[[rr_i]],
                                cft = cft[[cft_i]], weights = weight_matrix[, boot_i],
                                is_paf = is_paf)

      }
      #Transform to data.frame
      pif_replicates <- as.data.frame(pif_replicates)
      pif_replicates[,"counterfactual"] <- names(cft)[cft_i]
      pif_replicates[,"relative_risk"]  <- names(rr)[rr_i]

      return(pif_replicates)
    }
  }
  cli::cli_progress_done()


  #Create a pif_class object
  colnames(pif_values)[1] <- ifelse(is_paf,
                                    "population_attributable_fraction", "potential_impact_fraction")

  pif_obj <- pif_class(n_boot_simulations = ncol(weight_matrix),
                     is_paf = is_paf,
                     bootstrap_design = design,
                     theta_simulations = sim_theta,
                     pif_classulations = pif_values)

  return(pif_obj)
}


#' @title Pif from `data.frame`
#' @param df A data.frame in which to run the `pif` function.
#' @inheritParams pif
#' @keywords internal
#' @note In previous version of the package this function was called `pif.empirical`.
#' @seealso [pif()] [paf()]
#' @return A `vector` object containing the simulation for the `pif` (respectively `paf`),
#' the populational average relative risk and the populational average counterfactual.
internal_pif_data_frame <- function(df, theta, rr, cft, weights = NULL, is_paf = FALSE) {

  # Estimate weighted sums
  if (is.null(weights)) {
    mean_rr <- mean(as.matrix(rr(df, theta)))
  } else {
    mean_rr <- weighted.mean(as.matrix(rr(df, theta)), weights)
  }

  if (is_paf) {
    mean_cft <- 1
  } else {
    if (is.null(weights)) {
      mean_cft <- mean(as.matrix(rr(cft(df), theta)))
    } else {
      mean_cft <- weighted.mean(as.matrix(rr(cft(df), theta)), weights)
    }
  }

  # Calculate PIF
  if (is.infinite(mean_rr)) {
    cli::cli({
      cli::cli_alert_danger("Expected value of Relative Risk is not finite")
      cli::cli_rule("Tip:")
      cli::cli_text(
        "To debug you should call the relative risk function {.code rr} over your data.
        The problem is either:"
      )
      cli::cli_ul()
      cli::cli_li("The data {.code X} generates an infinite relative risk.
                  Try {.code rr(X, theta)} with the average parameters {.code theta}")
      cli::cli_li("The simulations for {.code theta} result in an infinite relative risk.
                  Try simulating {.code n_sim} thetas and then check for failures in:
                  {.code for (i in 1:n_sim) rr(X, theta[i])}")
      cli::cli_end()
      cli::cli_rule()
    })
  }

  if (is.infinite(mean_cft)) {
    cli::cli({
      cli::cli_alert_danger("Expected value of Relative Risk under counterfactual is not finite")
      cli::cli_rule("Tip:")
      cli::cli_text(
        "To debug you should call the relative risk function {.code rr} with the
        counterfactual {.code cft} over your data. The problem is either:"
      )
      cli::cli_ul()
      cli::cli_li("The counterfactual {.code cft} generates an infinite relative risk under the
                  counterfactual scenario. Try {.code rr(cft(X), theta)} with
                  the average {.code theta}")
      cli::cli_li("The simulations for {.code cft} result in an infinite relative risk under the
                  counterfactual scenario. Try simulating {.code n_sim} thetas and then check for
                  failures in: {.code for (i in 1:n_sim) rr(cft(X), theta[i])}")
      cli::cli_end()
      cli::cli_rule()
    })
  }

  pif <- 1 - mean_cft / mean_rr

  return(c("fraction" = pif, "average_relative_risk" = mean_rr,
           "average_counterfactual" = mean_cft))
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
#' @param rr (`function`/`list`) a relative risk function with two parameters: a `data.frame` called
#' `X` containing the individual-level exposure and covariates, and `theta` (in that order). It can
#' also be a list of several relative risk functions to apply with each function being a
#' different modelling scenario.
#'
#' @param cft (`function`/`list`) a counterfactual function that takes a `data.frame`, `X`, of
#' individual-level exposure and covariates and returns a new `data.frame` of individual-level
#' counterfactual exposure and covariates. It can also be a list of several counterfactual
#' functions to apply with each function being a different modelling scenario.
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
#' @param is_paf (`boolean`) Whether the function being estimated is the Population
#' Attributable Fraction (`is_paf = TRUE`) or the Potential Impact Fraction
#' (`is_paf = FALSE`)
#'
#' @param weights  (`vector`) If you are not following the recommended version and use a
#' `svydesign` object for the design you can still use `weights` to associate
#' weights to your estimation. Beware that it might not give accurate estimations of the variance
#' nor the uncertainty intervals.
#'
#' @param .options.future List of additional options for [doFuture::%dofuture%()].
#'
#' @param ... Additional parameters for [svrep::as_bootstrap_design()].
#'
#' @return A [pif_class()] object containing the bootstrap simulations for the
#' potential impact fraction, the average relative risk, and the average
#' counterfactual if applicable.
#'
#' @section Additional parallelization options:
#' Faster computation occurs when doing parallelization which allows to use more cores in your
#' machine. Parallelization utilizes  the [future::future()] package. For paralelization to work you
#' need to establish a plan (see [future::plan()] for more information). The most common
#' way to create parallelization in your local machine is to do:
#' ```
#' plan(multisession) #Start parallelization
#' pif(...)
#' plan(sequential) #Return to 'normal'
#' ```
#'
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
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10,
#' )
#'
#' # EXAMPLE 2
#' # Now do the same but using a replicate design
#' rep_design <- svrep::as_bootstrap_design(design, replicates = 10)
#' pif(rep_design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03)
#' )
#'
#' # EXAMPLE 3
#' # Calculate two different relative risks
#' rr <- list(
#'  function(X, theta) {exp(theta[1] * X[, "systolic_blood_pressure"] / 100)},
#'  function(X, theta) {exp(theta[2] * X[, "systolic_blood_pressure"] / 100 + theta[3]* X[,"age"])}
#' )
#'
#' # Calculate three counterfactual scenarios of SBP reduction from 1 to 3 mmhg
#' cft <- list(
#' function(X){X[, "systolic_blood_pressure"]  <- X[, "systolic_blood_pressure"]  - 1; return(X)},
#' function(X){X[, "systolic_blood_pressure"]  <- X[, "systolic_blood_pressure"]  - 2; return(X)},
#' function(X){X[, "systolic_blood_pressure"]  <- X[, "systolic_blood_pressure"]  - 3; return(X)}
#' )
#'
#' pif(design,
#'   theta = log(c(1.05, 1.38, 1.21)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03, 0.025),
#'   n_bootstrap_samples = 10,
#' )
#' @references \insertAllCited{}
#' @seealso [paf()] [plot.pif_class()] [summary.pif_class()]
#' @keywords pif impact-fraction paf
#' @export
#' @importFrom doFuture  %dofuture%
#' @md

pif <- function(design,
                theta,
                rr,
                cft = NULL,
                additional_theta_arguments,
                n_bootstrap_samples = NULL,
                theta_distribution = "default",
                is_paf = FALSE,
                weights = NULL,
                .options.future = list(seed = TRUE),
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
                         is_paf = is_paf, .options.future = .options.future,...)
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
#' @return A [pif_class()] object containing the bootstrap simulations for the
#' population attributable fraction, and the average relative risk.
#'
#' @inheritSection pif Additional parallelization options
#'
#' @examples
#' # Use the ensanut dataset
#' data(ensanut)
#'
#' # EXAMPLE 1
#' # Setup the survey design
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#' rr <- function(X, theta) {
#'   exp(
#'     theta[1] * X[, "age"] + theta[2] * X[, "systolic_blood_pressure"] / 100)
#' }
#' paf(design,
#'   theta = log(c(1.05, 1.38)), rr,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10,
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
#' @seealso [pif()] [plot.pif_class()]
#' @keywords paf attributable-fraction
#' @export
#' @importFrom doFuture  %dofuture%
#' @md

paf <- function(design,
                theta,
                rr,
                additional_theta_arguments,
                n_bootstrap_samples = NULL,
                theta_distribution = "default",
                weights = NULL,
                .options.future = list(seed = TRUE),
                ...) {


  return(
    pif(design = design,
        theta = theta,
        rr = rr,
        additional_theta_arguments = additional_theta_arguments,
        n_bootstrap_samples = n_bootstrap_samples,
        theta_distribution = theta_distribution,
        is_paf = TRUE,
        .options.future = .options.future,
        weights = weights,
         ...)
  )
}
