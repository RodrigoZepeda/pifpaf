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
#' @param return_pif_replicates (`boolean`) Whether to return the simulated impact fractions
#' or not.
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
#'\mjdeqn{
#' \Big[\widehat{\text{PIF}}_{\text{B},\alpha/2}, \widehat{\text{PIF}}_{\text{B},1-\alpha/2}\Big]
#'}{quantile(pif, c(alpha/2, 1-alpha/2))}
#' where \mjseqn{\widehat{\text{PIF}}_{\text{B},k}} represents the kth sample quantile
#' from the bootstrap simulation of the potential impact fraction estimator.
#'
#' @section Additional parallelization options:
#' By default the function uses [foreach::foreach] to parallelize creating a cluster
#' via [doParallel::registerDoParallel]. For finner parallelization control you can do:
#' ```
#' #Create the cluster outside using whatever you want
#' myCluster <- makeCluster(3, type = "PSOCK")
#'
#' #Register the cluster
#' registerDoParallel(myCluster)
#'
#' #Call the function using parallel = TRUE to use %dopar% and num_cores = NULL to avoid setting
#' pif_survey_bootstrap(..., parallel = TRUE, num_cores = NULL)
#'
#' #Stop the cluster
#' stopCluster(myCluster)
#' ```
#' @examples
#' #Use the ensanut dataset
#' data(ensanut)
#'
#' #EXAMPLE 1
#' #Setup the survey design
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#' rr <- function(X, theta){exp(-2 + theta[1]*X[,"age"] + theta[2]*X[,"systolic_blood_pressure"]/100)}
#' cft <- function(X){X[,"systolic_blood_pressure"] <- X[,"systolic_blood_pressure"] - 5; return(X)}
#' pif_survey_bootstrap(design, theta = log(c(1.05, 1.38)), rr, cft, additional_theta_arguments = c(0.01, 0.03),
#'  n_bootstrap_samples = 10, parallel = F,
#' )
#'
#' #EXAMPLE 2
#' #Now do the same but using a replicate design
#' options(survey.lonely.psu = "adjust")
#' rep_design <- svrep::as_bootstrap_design(design, replicates = 10)
#'
#' @references \insertAllCited{}
#'
#' @seealso [paf()], [pif()], [pif_approximate()], [pif_delta_method()].
#'
#' @keywords internal, pif, impact-fraction
#' @export
#' @importFrom parallel detectCores
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @md

pif_survey_bootstrap <- function(design,
                                 theta,
                                 rr,
                                 cft,
                                 additional_theta_arguments,
                                 n_bootstrap_samples = NULL,
                                 theta_distribution = "default",
                                 uncertainty_interval_type = c("wald","percentile","none"),
                                 parallel = TRUE,
                                 num_cores = 1,
                                 confidence_level = 0.95,
                                 return_replicates = FALSE,
                                 is_paf = FALSE,
                                 ...){

  #Validators:
  .design <- validate_survey_design(design = design, n_bootstrap_samples = n_bootstrap_samples, ...)
  .confidence_level <- validate_confidence_level(confidence_level = confidence_level)
  .num_cores <- validate_number_of_cores(num_cores = num_cores)

  #Parallel setup
  `%docommand%` <-ifelse(parallel & .num_cores > 1, foreach::`%dopar%`, foreach::`%do%`)

  if(!is.null(.num_cores) & .num_cores > 1){
    cl <- parallel::makeCluster(.num_cores)
    doParallel::registerDoParallel(cl)
  }

  #Get the weights
  #https://stackoverflow.com/questions/73627746/how-do-i-get-the-weights-from-a-survey-design-object-in-r
  .data <- .design[["variables"]]
  .weight_matrix <- stats::weights(.design, type = "analysis")

  #Get the simulated thetas
  if (!is.function(theta_distribution) && theta_distribution == "default" ){
    ifelse(
      (length(additional_theta_arguments) == 1),

      #Sigma as a variable
      .theta_args <- list("n" = ncol(.weight_matrix),
                          "sigma" = as.matrix(additional_theta_arguments[[1]]),
                          "mean" = theta),
    ifelse(
      is.vector(additional_theta_arguments),

      #Diagonal matrix for independent
      .theta_args <- list("n" = ncol(.weight_matrix),
                          "sigma" = diag(additional_theta_arguments),
                          "mean" = theta),

      #Default (matrix)
      .theta_args <- list("n" = ncol(.weight_matrix),
                          "sigma" = additional_theta_arguments,
                          "mean" = theta)

    ))

    .theta_distribution <- mvtnorm::rmvnorm

  } else {
    #Append
    .theta_args <- append(.theta_args, list("n" = ncol(.weight_matrix)))
    .theta_distribution <- theta_distribution
  }

  #Call function
  .sim_theta <- do.call(.theta_distribution, args = .theta_args)

  #Loop through the foreach
  .pif_replicates <- foreach::foreach(sim = 1:ncol(.weight_matrix), .combine = c, .inorder = FALSE,
                                      .multicombine	= TRUE,
                                      .noexport = c("design")) %docommand% {

    #Obtain the pif
    pif_data_frame(df = .data, theta = .sim_theta[sim,], rr = rr, cft = cft,
                   weights = .weight_matrix[,sim], is_paf = is_paf)
  }

  if(!is.null(.num_cores) & .num_cores > 1){
    stopImplicitCluster(cl)
  }

  #Get point estimate
  .point_estimate <- mean(.pif_replicates)
  .std_pif <- sd(.pif_replicates)

  #Get alpha of Uncertainty interval
  .alpha_ui <- 1 - .confidence_level

  if (uncertainty_interval_type[1] == "percentile"){
    .interval <- quantile(.pif_replicates, c(.alpha_ui/2, 1 - .alpha_ui/2))
    names(.interval) <- c("Lower", "Upper")
  } else if (uncertainty_interval_type[1] == "wald"){
    .interval <- c(
      "Lower" = .point_estimate - qt(1 - .alpha_ui/2, df = ncol(.weight_matrix))*.std_pif,
      "Upper" = .point_estimate + qt(1 - .alpha_ui/2, df = ncol(.weight_matrix))*.std_pif
    )
  }

  .pif_simulations <- list("Point" = .point_estimate,
               "Interval" = .interval,
               "Confidence" = .confidence_level,
               "Standard_Deviation" = .std_pif)

  if (return_replicates){
    .pif_simulations <- append(.pif_simulations, list("Replicates" = .pif_replicates))
  }

  return(.pif_simulations)

}


#' @title Pif from dataframe
#'
#' @note In previous version of the package this function was called `pif.empirical`.
pif_data_frame <- function(df, theta, rr, cft,
                           weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                           is_paf = FALSE){


  #Estimate weighted sums
  .mux   <- weighted.mean(as.matrix(rr(df, theta)), weights)
  if (is_paf){
    .mucft <- 1
  } else {
    .mucft <- weighted.mean(as.matrix(rr(cft(df), theta)), weights)
  }

  #Calculate PIF
  if( is.infinite(.mux) ){
    cli::cli_alert_warning("Expected value of Relative Risk is not finite")
  }

  if( is.infinite(.mucft) ){
    cli::cli_alert_warning("Expected value of Relative Risk under counterfactual is not finite")
  }

  .pif   <- 1 - .mucft/.mux

  return(.pif)

}


