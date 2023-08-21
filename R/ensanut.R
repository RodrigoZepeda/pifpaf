#' @title ENSANUT 2018
#'
#' @description
#' Data on the Mexican [National Health and Nutrition Survey of 2018 (ENSANUT-2018)](https://ensanut.insp.mx/).
#'
#' @format Information contains: A tibble containing the following information
#' * `age` Age of the individual (years).
#' * `sex` Sex of the individual (`Male` or `Female`).
#' * `weight` Survey weights associated to the survey.
#' * `strata` Strata of the survey.
#' * `hypertension` If individual has hypertension `TRUE`, `FALSE` otherwise.
#' * `delta_na_phase_1` Proposed change in sodium consumption (counterfactual).
#' * `systolic_blood_pressure` Mean systolic blood pressure of the individual registered by survey.
#' * `age_group` Age group of the individual
#'
#' @examples
#' ensanut
#' #Create a survey design as follows:
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
"ensanut"
