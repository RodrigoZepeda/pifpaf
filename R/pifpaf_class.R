#' @title pif_sim class
#'
#' @description An [S7::S7_class()] for saving the bootstrap simulations of the potential
#' impact fraction [pif()] and the population attributable fraction [paf()].
#'
#' @note The use of this class are for development only.
#'
#' @param pif_simulations Bootstrap simulations for the impact fraction or the attributable
#' fraction.
#' @param bootstrap_design The survey bootstrap design utilized for estimating `pif` or `paf`.
#' @param theta_simulations Simulations for the parameter `theta` of the relative risk.
#' @inheritParams pif
#' @export
pif_sim <- S7::new_class("pif_sim",
                         properties = list(
                           is_paf            = S7::class_logical,
                           bootstrap_design  = S7::class_any,
                           theta_simulations = S7::class_any,
                           pif_simulations   = S7::class_data.frame
                         ))


#' @title Extract point estimate of potential impact fraction or population attributable fraction
#' @description Method for obtaining the point estimate of a `pif_sim` object (i.e. the
#' point-estimate of the potential impact fraction or population attributable fraction via
#' the bootstrap).
#' @param x A `pif_sim` object
#' @param ... additional arguments for [base::mean()]
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' coef(pifsim)
#' @export
coef.pif_sim <- function(x, ...) {
  vars  <- c("counterfactual", "relative_risk")
  dplyr::summarise(S7::prop(x,"pif_simulations"),
                   dplyr::across(.cols = -dplyr::any_of(vars),
                                 .fns = ~mean(., ...)),
                                 .by = vars)
}


#' @title Extract confidence intervals for the potential impact fraction or the
#' population attributable fraction
#'
#' @description Method for obtaining confidence intervals for a `pif_sim` object
#'
#' @param x A `pif_sim` object
#'
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#'
#' @param method (`string`) either `'wald'` (recommended) or `'percentile'` for
#' uncertainty intervals based upon the bootstrap's percentile or the Wald-type approximation
#' (see confidence interval section)
#'
#' @param level (`double`) Confidence level for the uncertainty interval.
#' Defaults to `0.95`.
#'
#' @param ... additional methods to pass to [mean()] and [stats::sd()]
#' if `wald` or [stats::quantile()] if `percentile` type intervals
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
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' confint(pifsim)
#' @export
confint.pif_sim <- function(x, parm = NULL, level = 0.95, method = c("wald", "percentile"),
                            ...) {

  #Validate the confidence level
  confidence_level    <- validate_confidence_level(confidence_level = level)

  # Get alpha of uncertainty interval
  alpha_ui <- 1 - confidence_level

  #Number of simulations
  sims <- S7::prop(x,"pif_simulations")
  if (!is.null(parm)){
    sims <- sims[,c(parm, "counterfactual", "relative_risk")]
  }

  #Chech the method
  vars  <- c("counterfactual", "relative_risk")

  #Get the quantiles
  prob  <- c(alpha_ui/2, 1 - alpha_ui/2)

  if (method[1] == "percentile") {

    interval <- cbind(
      dplyr::reframe(sims,
                     dplyr::across(.cols = -dplyr::any_of(vars),
                                   .fns = ~stats::quantile(., prob = !!prob, ...)),
                     .by = vars),
      "type" = paste(c("Lower", "Upper"), paste0(round(prob*100, 2), "%"))
    )

  } else if (method[1] == "wald") {

    #Get the size of the bootstrap
    nvals <- ncol(stats::weights(S7::prop(x,"bootstrap_design")))

    # Get point estimate
    cifun <- function(x, ...){
      c(
        "Lower" = mean(x, ...) - stats::qt(1 - alpha_ui / 2, df = nvals)*sd(x, ...),
        "Upper" = mean(x, ...) + stats::qt(1 - alpha_ui / 2, df = nvals)*sd(x, ...)
      )
    }

    interval <- cbind(
      dplyr::reframe(sims,
                     dplyr::across(.cols = -dplyr::any_of(vars),
                                   .fns = !!cifun), .by = vars),
      "type" = paste(c("Lower", "Upper"), paste0(round(prob*100, 2), "%"))
    )

  } else {
    cli::cli_abort("Invalid `method` for confidence intervals.
                   Please choose `'wald'` or `'percentile'`")
  }

  return(interval)
}

#' @title Extract mean of bootstrap simulations of potential impact fraction or population
#' attributable fraction
#' @description Method for obtaining the arithmetic mean of bootstrap simulations
#' of a `pif_sim` object (i.e. the point-estimate of the potential impact fraction or
#' population attributable fraction via the bootstrap).
#' @param x A `pif_sim` object
#' @param ... additional arguments to pass to [base::mean()]
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' mean(pifsim)
#'
#' @export
mean.pif_sim <- function(x, ...) {
  vars  <- c("counterfactual", "relative_risk")
  dplyr::summarise(S7::prop(x,"pif_simulations"),
                   dplyr::across(.cols = -dplyr::any_of(vars),
                                 .fns = ~mean(., ...)),
                   .by = vars)
}

#' @title Extract standard deviations of bootstrap simulations of potential impact fraction
#' or population attributable fraction
#' @description Method for obtaining the standard deviartion of bootstrap simulations
#' of a `x` object (i.e. the point-estimate of the potential impact fraction or
#' population attributable fraction via the bootstrap).
#' @param x A `pif_sim` object
#' @param ... additional arguments to pass to [stats::sd()]
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' sd_bootstrap(pifsim)
#' @export
sd_bootstrap <- function(x, ...) {
  vars  <- c("counterfactual", "relative_risk")
  dplyr::summarise(S7::prop(x,"pif_simulations"),
                   dplyr::across(.cols = -dplyr::any_of(vars),
                                 .fns = ~stats::sd(., ...)),
                   .by = vars)
}

#' @title Extract the variance of bootstrap simulations of potential impact fraction
#' or population attributable fraction
#' @description Method for obtaining the variance of bootstrap simulations
#' of a `x` object (i.e. the point-estimate of the potential impact fraction or
#' population attributable fraction via the bootstrap).
#' @param x A `pif_sim` object
#' @param ... additional arguments to pass to [stats::var()]
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' var_bootstrap(pifsim)
#' @export
var_bootstrap <- function(x, ...) {
  vars  <- c("counterfactual", "relative_risk")
  dplyr::summarise(S7::prop(x,"pif_simulations"),
                   dplyr::across(.cols = -dplyr::any_of(vars),
                                 .fns = ~stats::var(., ...)),
                   .by = vars)
}

#' @title Extract quantiles of the potential impact fraction or
#' population attributable fraction bootstrap simulations
#' @description Method for obtaining the variance-covariance matrix of bootstrap simulations
#' of a `pif_sim` object (i.e. the point-estimate of the potential impact fraction or
#' population attributable fraction via the bootstrap).
#' @param x A `pif_sim` object
#' @param prob numeric vector of probabilities with values between 0 and 1
#' @param ... additional arguments to pass to [stats::quantile]
#' @export
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' quantile_bootstrap(pifsim, c(0.1, 0.3, 0.5))
#'
quantile_bootstrap <- function(x, prob, ...) {
  vars  <- c("counterfactual", "relative_risk")
  cbind(
    dplyr::reframe(S7::prop(x,"pif_simulations"),
                   dplyr::across(.cols = -dplyr::any_of(vars),
                                 .fns = ~stats::quantile(., prob = !!prob, ...)),
                     .by = vars),
    "quantile" = prob
  )
}

#' @title Extract number of bootstrap simulations
#' @description Method for obtaining the number of bootstrap simulations used to construct
#' a `pif_sim` object.
#' @param x A `pif_sim` object
#' @export
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' n_bootstrap(pifsim)
#'
n_bootstrap <- function(x) {
  nrow(S7::prop(x,"pif_simulations"))
}

#' @title Return potential impact fraction or population attributable fraction summary
#' @description Method for obtaining the summary of bootstrap simulations
#' of a `pif_sim` object (i.e. the point-estimate of the potential impact fraction or
#' population attributable fraction via the bootstrap).
#' @param object A `pif_sim` object
#' @inheritParams confint.pif_sim
#' @param ... Additional methods to pass to [coef.pif_sim()]
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' summary(pifsim)
#' @export
summary.pif_sim <- function(object, parm = NULL, level = 0.95, method = c("wald", "percentile"),
                            ...) {
  list("Type" = get_fraction_type(object),
    "Number_of_bootstrap_simulations" = n_bootstrap(object),
    "Point_estimates" = coef.pif_sim(object, ...),
    "Confidence_intervals" = confint.pif_sim(object,
                                             parm = parm, level = level, method = method)
  )
}

#' @title Get whether object is population attributable fraction or potential impact fraction
#' @description Returns whether a `pif_sim` object was specified by the creator as population
#' attributable fraction or potential impact fraction
#' @param object A `pif_sim` object
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' get_fraction_type(pifsim)
#' @export
get_fraction_type <- function(object){
  ifelse(S7::prop(object,"is_paf"),
                "Population Attributable Fraction (PAF)", "Potential Impact Fraction (PIF)")
}

#' @title Get bootstrap simulations from `pif` and `paf`
#' @description Returns the bootstrap simulations for `pif` and `paf`
#' @param object A `pif_sim` object
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' get_bootstrap_simulations(pifsim)
#' @export
get_bootstrap_simulations <- function(object){
  S7::prop(object,"pif_simulations")
}

#' @title Get bootstrap simulations for `theta`
#' @description Returns the bootstrap simulations for the parameter `theta` of `pif` or `paf`
#' @param object A `pif_sim` object
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' get_theta_simulations(pifsim)
#' @export
get_theta_simulations <- function(object){
  S7::prop(object,"theta_simulations")
}

#' @title Transform to data.frame
#' @description Method for transforming a potential impact fraction or a
#' population attributable fraction to a data.frame
#' @param x A `pif_sim` object
#' @param ...  Additional methods to pass to [summary.pif_sim()]
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' as.data.frame(pifsim)
#' @export
as.data.frame.pif_sim <- function(x, ...){
  summary_data <- summary.pif_sim(x, ...)
  summary_data[["Point_estimates"]][,"type"] <- "point_estimate"

  rbind(summary_data[["Point_estimates"]], summary_data[["Confidence_intervals"]])
}

#' @title Print an impact or attributable fraction
#' @description Method for printing a potential impact fraction or a population attributable fraction
#' @param x A `pif_sim` object
#' @param max_sim_print Maximum number of simulations to print
#' @param ... Additional arguments to pass to [summary.pif_sim]
#' @examples
#' #Example 1
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
#' pifsim <- pif(design,
#'   theta = log(c(1.05, 1.38)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10, parallel = FALSE
#' )
#' print(pifsim)
#' print(pifsim, max_sim_print = 1)
#' @export
print.pif_sim <- function(x, max_sim_print = 10000, ...) {
  cli::cli_rule(get_fraction_type(x))
  if (nrow(S7::prop(x, "pif_simulations")) <= max_sim_print){
    print(dplyr::select(as.data.frame(x),
                        c(dplyr::starts_with("potential"),
                          dplyr::starts_with("population"),
                          dplyr::starts_with("type"),
                          dplyr::starts_with("average"),
                          dplyr::starts_with("relative"),
                          dplyr::starts_with("counterfactual"))), ...)
  } else {
    cli::cli_alert_info("Too many simulations to print.")
  }

  cli::cli({
    cli::cli_rule()
    cli::cli_ul()
    cli::cli_li("Number of bootstrap simulations: {.value {n_bootstrap(x)}}")
    if (n_bootstrap(x) < 1000){
      cli::cli_alert_danger(
        "A low number of bootstrap simulations will result in an unstable estimate."
      )
    }
    cli::cli_li("Use {.code as.data.frame} to access values.")
    cli::cli_li("Use {.code summary} to save list of main results.")
    cli::cli_end()
  })
}
