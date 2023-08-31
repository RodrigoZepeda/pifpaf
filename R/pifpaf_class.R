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
                           rr                        = S7::class_function,
                           cft                       = S7::class_function,
                           is_paf                    = S7::class_logical,
                           bootstrap_design          = S7::class_any,
                           theta_simulations         = S7::class_any,
                           pif_simulations           = S7::class_numeric,
                           uncertainty_interval_type = S7::class_character,
                           confidence_level          = S7::class_numeric
                         ))


#' @title Extract point estimate of potential impact fraction or population attributable fraction
#' @description Method for obtaining the point estimate of a `pif_sim` object (i.e. the
#' point-estimate of the potential impact fraction or population attributable fraction via
#' the bootstrap).
#' @param pif_sim A `pif_sim` object
#' @param ... additional arguments
#' @export
#' @examples
#' some_disease <- pif_sim(pif_simulations = rbeta(100, 0.1, 0.2), is_paf = FALSE)
#' coef(some_disease)
#'
coef.pif_sim <- S7::new_generic("coef.pif_sim", dispatch_args = "pif_sim")

S7::method(coef.pif_sim, pif_sim) <- function(pif_sim, ...) {
  ifelse(pif_sim@is_paf,
         c("paf" = mean(pif_sim@pif_simulations)),
         c("pif" = mean(pif_sim@pif_simulations)))
}

#' @title Extract confidence intervals for the potential impact fraction or the
#' population attributable fraction
#' @description Method for obtaining confidence intervals for a `pif_sim` object
#' @param pif_sim A `pif_sim` object
#' @param ... additional arguments
#' @export
#' @examples
#' some_disease <- pif_sim(pif_simulations = rbeta(100, 0.1, 0.2), is_paf = FALSE)
#' coef(some_disease)
#'
confint.pif_sim <- S7::new_generic("confint.pif_sim", dispatch_args = "pif_sim")

S7::method(confint.pif_sim, pif_sim) <- function(pif_sim, parm = NULL, level = 0.95, ...) {
  print("Confident")
}

#TODO: print, mean, quantiles, variance, vcov, summary
#TODO: save the rr and the cft simulations separately
#TODO: create a multiple pif object
