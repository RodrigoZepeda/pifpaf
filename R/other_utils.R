#' Mean estimator's standard error from asymptotically normal confidence interval
#'
#' @description Auxiliary function to estimate the standard error (variance divided by sample size)
#' from an asymptotically normal confidence interval of level `1 - alpha` by implementing
#' the following function:
#'
#' \loadmathjax
#' \mjdeqn{
#'  \text{se} = \dfrac{U - L}{2\cdot Z_{1 - \alpha/2}}
#' }{se = (U - L)/2*qnorm(1 - alpha/2)}
#' where `U` and `L` are the upper and lower bounds of the interval
#'
#' @param upper_limit Upper bound for the confidence interval
#' @param lower_limit Lower bound for the confidence interval
#' @param confidence_level Confidence level of the confidence interval
#'
#' @examples
#' # Calculate the confidence interval for a normal mean
#' x <- rnorm(100)
#' interval_mu_x <- t.test(x)$conf.int
#' estimator_se_from_ci(interval_mu_x[1], interval_mu_x[2], 0.95)
#'
#' @returns The standard error of the estimator
#'
#' @seealso [estimator_var_from_ci()]
#' @export
estimator_se_from_ci <- function(upper_limit, lower_limit, confidence_level = 0.95){

  #Get correct confidence level
  confidence_level <- ifelse(confidence_level > 1, confidence_level/100, confidence_level)

  if(confidence_level > 1 | confidence_level < 0){
    cli::cli_abort(
      "Invalid confidence interval level not between `0` and `1`. Maybe a typo?"
    )
  }

  # Calculate the Z-score based on the confidence level
  z_score <- stats::qnorm((1 + confidence_level) / 2)

  # Calculate the standard error
  standard_error <- (upper_limit - lower_limit) / (2 * z_score)

  return(standard_error)
}

#' Mean estimator's variance from asymptotically normal confidence interval
#'
#' @description Auxiliary function to estimate the estimator's variance (_i.e. standard error
#' squared) from an asymptotically normal confidence interval of level `1 - alpha` by implementing
#' the following function:
#'
#' \loadmathjax
#' \mjdeqn{
#'  \text{var} = \Bigg(\dfrac{U - L}{2\cdot Z_{1 - \alpha/2}}\Bigg)^2
#' }{se = (U - L)/2*qnorm(1 - alpha/2)}
#' where `U` and `L` are the upper and lower bounds of the interval
#'
#' @param upper_limit Upper bound for the confidence interval
#' @param lower_limit Lower bound for the confidence interval
#' @param confidence_level Confidence level of the confidence interval
#'
#' @examples
#' # Calculate the confidence interval for a normal mean
#' x <- rnorm(100)
#' interval_mu_x <- t.test(x)$conf.int
#' estimator_var_from_ci(interval_mu_x[1], interval_mu_x[2], 0.95)
#'
#' @returns The variance of the estimator
#' @note Remember that the variance of the mean-estimator from asymptotically normal
#' confidence intervals is the estimator's standard error squared
#'
#' @seealso [estimator_se_from_ci()]
#' @export
estimator_var_from_ci <- function(upper_limit, lower_limit, confidence_level = 0.95){
  estimator_se_from_ci(upper_limit = upper_limit, lower_limit = lower_limit,
                       confidence_level = confidence_level)^2
}
