#' @title Joint PIF/PAF type
#'
#' @description
#' Function that determines whether a collection of potential impact fractions or
#' population attributable fractions are pif or paf
#'
#' @param ... `pif_class` objects (_i.e._ estimations of `pif` or `paf`)
#'
#' @references \insertAllCited{}
#' @examples
#' # Use the ensanut dataset
#' data(ensanut)
#'
#' # EXAMPLE 1
#' # Setup the survey design
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#' rr_1 <- function(X, theta) {
#'   exp(-2 +
#'     theta[1] * X[, "age"] + theta[2] * X[, "systolic_blood_pressure"] / 100)
#' }
#'
#' rr_2 <- function(X, theta) {
#'   exp(-1 +
#'     theta[1] * X[, "age"] )
#' }
#'
#' rr_3 <- function(X, theta) {
#'   exp(-3 +
#'     theta[1] * X[, "age"] )
#' }
#'
#' cft <- function(X) {
#'   X[, "systolic_blood_pressure"] <- X[, "systolic_blood_pressure"] - 5
#'   return(X)
#' }
#'
#' pif_1 <- pif(design,
#'   theta = log(c(1.25, 1.68)), rr_1, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10,
#' )
#'
#' pif_2 <- pif(design,
#'   theta = log(c(1.12, 1.45)), rr_2, cft,
#'   additional_theta_arguments = c(0.02, 0.025), n_bootstrap_samples = 10,
#' )
#'
#' pif_3 <- pif(design,
#'   theta = log(c(2.25, 2.57)), rr_3, cft,
#'   additional_theta_arguments = c(0.01, 0.025), n_bootstrap_samples = 10,
#' )
#'
#' overall_fraction_type(pif_1, pif_2, pif_3)
#'
#' @seealso [pif()] [paf()]
#' @return A `pif_aggregated` object estimating the potential impact fraction for
#' individual-level data from multiple distinct impact fractions.
#' @export

overall_fraction_type <- function(...){
  #Get all pif values
  pifs <- list(...)
  if (do.call(validate_pifs, pifs)){

    if (any(sapply(pifs, get_fraction_type) == "Potential Impact Fraction (PIF)")){
      type_pif <- "Potential Impact Fraction (PIF)"
    } else {
      type_pif <- "Population Attributable Fraction (PAF)"
    }

  } else {
    cli::cli_abort(
      "Argument is not `pif_class`. Did you estimate it using the `pif`/`paf` functions?"
    )
  }
  return(type_pif)
}

#' @title Joint PIF: Combine impact fractions from independent and uncorrelated risk factors
#'
#' @description
#' A function to combine multiple  impact fractions from independent uncorrelated risk factors. The
#' function follows the formula from \insertCite{ezzati2003estimates;textual}{pifpaf}:
#' \loadmathjax
#' \mjdeqn{
#'  \text{PIF} = 1 - \prod\limits_{k} \big(1 - \text{PIF}_i)
#' }{pif = 1 - (1 - pif[1])*(1 - pif[2])*...*(1 - pif[n])}
#'
#' @param pif_class_1 A `pif_class` object (_i.e._ an estimation of `pif` or `paf`)
#' @param pif_class_2 Another `pif_class` object (_i.e._ an estimation of a different `pif` or `paf`)
#' @param ... Additional `pif_class` objects (_i.e._ estimations of `pif` or `paf`)
#'
#' @references \insertAllCited{}
#' @examples
#' # Use the ensanut dataset
#' data(ensanut)
#'
#' # EXAMPLE 1
#' # Setup the survey design
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#' rr_1 <- function(X, theta) {
#'   exp(-2 +
#'     theta[1] * X[, "age"] + theta[2] * X[, "systolic_blood_pressure"] / 100)
#' }
#'
#' rr_2 <- function(X, theta) {
#'   exp(-1 +
#'     theta[1] * X[, "age"] + X[, "systolic_blood_pressure"] )
#' }
#'
#' rr_3 <- function(X, theta) {
#'   exp(-1 +
#'     theta[1] *X[, "systolic_blood_pressure"] )
#' }
#'
#' cft <- function(X) {
#'   X[, "systolic_blood_pressure"] <- X[, "systolic_blood_pressure"] - 5
#'   return(X)
#' }
#'
#' pif_1 <- pif(design,
#'   theta = log(c(1.25, 1.68)), rr_1, cft,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10,
#' )
#'
#' pif_2 <- pif(design,
#'   theta = log(c(1.12, 1.45)), rr_2, cft,
#'   additional_theta_arguments = c(0.02, 0.025), n_bootstrap_samples = 10,
#' )
#'
#' pif_3 <- pif(design,
#'   theta = 0.24, rr_3, cft,
#'   additional_theta_arguments = 0.01, n_bootstrap_samples = 10,
#' )
#'
#' pif_combine(pif_1, pif_2, pif_3)
#'
#' @seealso [pif()] [paf()]
#' @return A `pif_aggregated` object estimating the potential impact fraction for
#' individual-level data from multiple distinct impact fractions.
#' @export

pif_combine <- function(pif_class_1, pif_class_2, ...){
  #Get all pif values
  pifs <- list(pif_class_1, pif_class_2, ...)
  if (do.call(validate_pifs, pifs)){

    frac_name <- ifelse(do.call(overall_fraction_type, pifs) == "Potential Impact Fraction (PIF)",
                        "potential_impact_fraction", "population_attributable_fraction")

    #Create a list of simulations of all objects pif
    pifmat <- lapply(pifs, get_pif_simulations)
    i <- 0
    suppressMessages({
      pifmat <- foreach::foreach(i = 1:length(pifmat), .combine = dplyr::left_join) %do% {
          colnames(pifmat[[i]])[1] <- paste0(frac_name, "_", i)
          pifmat[[i]][,"sim"]      <- 1:nrow(pifmat[[i]])
          return(pifmat[[i]])
      }
    })

    #Now calculate the product of 1 - pif grouped by counterfactual and relative risk
    pifmat <- dplyr::transmute(
      dplyr::rowwise(
        dplyr::group_by(pifmat, !!as.symbol("counterfactual"), !!as.symbol("relative_risk"))),
            !!as.symbol(frac_name) :=  1 - prod(1 - dplyr::c_across(dplyr::starts_with(!!frac_name)))
    )

    pifmat <- dplyr::ungroup(pifmat)

    #Create a pif object
    pif_aggregated <- pif_class(
      n_boot_simulations = nrow(pifmat),
      is_paf = do.call(overall_fraction_type, pifs) == "Population Attributable Fraction (PAF)",
      bootstrap_design = NA,
      theta_simulations = NA_real_,
      pif_classulations = pifmat
    )

  } else {
    cli::cli_abort(
      "Argument is not `pif_class`. Did you estimate it using the `pif`/`paf` functions?"
    )
  }
  return(pif_aggregated)
}

#' @title Joint PAF: Combine attributable fractions from independent and uncorrelated risk factors
#'
#' @description
#' A function to combine multiple attributable fractions from independent uncorrelated risk factors. The
#' function follows the formula from \insertCite{ezzati2003estimates;textual}{pifpaf}:
#' \loadmathjax
#' \mjdeqn{
#'  \text{PAF} = 1 - \prod\limits_{k} \big(1 - \text{PAF}_i)
#' }{paf = 1 - (1 - paf[1])*(1 - paf[2])*...*(1 - paf[n])}
#'
#' @param paf_class_1 A `pif_class` object (_i.e._ an estimation of `paf`)
#' @param paf_class_2 Another `pif_class` object (_i.e._ an estimation of a different `paf`)
#' @param ... Additional `pif_class` objects (_i.e._ estimations of `paf`)
#'
#' @references \insertAllCited{}
#' @examples
#' # Use the ensanut dataset
#' data(ensanut)
#'
#' # EXAMPLE 1
#' # Setup the survey design
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#' rr_1 <- function(X, theta) {
#'   exp(-2 +
#'     theta[1] * X[, "age"] + theta[2] * X[, "systolic_blood_pressure"] / 100)
#' }
#'
#' rr_2 <- function(X, theta) {
#'   exp(-1 +
#'     theta[1] * X[, "age"] + X[, "systolic_blood_pressure"] )
#' }
#'
#' rr_3 <- function(X, theta) {
#'   exp(-1 +
#'     theta[1] *X[, "systolic_blood_pressure"] )
#' }
#'
#' paf_1 <- paf(design,
#'   theta = log(c(1.25, 1.68)), rr_1,
#'   additional_theta_arguments = c(0.01, 0.03), n_bootstrap_samples = 10,
#' )
#'
#' paf_2 <- paf(design,
#'   theta = log(c(1.12, 1.45)), rr_2,
#'   additional_theta_arguments = c(0.02, 0.025), n_bootstrap_samples = 10,
#' )
#'
#' paf_3 <- paf(design,
#'   theta = 0.24, rr_3,
#'   additional_theta_arguments = 0.01, n_bootstrap_samples = 10,
#' )
#'
#' paf_combine(paf_1, paf_2, paf_3)
#'
#' @seealso [pif()] [paf()]
#' @return A `paf_aggregated` object estimating the population attributable fraction for
#' individual-level data from multiple distinct attributable fractions.
#' @export
paf_combine <- function(paf_class_1, paf_class_2, ...){
  ftype <- do.call(overall_fraction_type, list(paf_class_1, paf_class_2, ...))
  if (ftype != "Population Attributable Fraction (PAF)"){
    cli::cli_abort(
      "Object is not a population attributable fraction (PAF). Use `pif_combine` instead.")
  } else {
    return(pif_combine(paf_class_1, paf_class_2, ...))
  }
}
