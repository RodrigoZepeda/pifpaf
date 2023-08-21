#' @title Validate confidence level
#'
#' @description A function that takes the a confidence level as input. It can be inputed as:
#' * A number between `0` and `1`.
#' * A number between `1` and `100` which is then divided by 100.
#' * Any other input generates an error.
#'
#' @param confidence_level (`double`) Confidence level for the uncertainty interval. If its between `0`
#' and `1` it uses it _as is_. If it is between `1` and `100` it divides it by `100` to
#' transform it to a probability. Defaults to `0.95`.
#'
#' @returns (`double`) A list containing a confidence level between `0` and `1`, and the values
#' of `alpha` and `1 - alpha/2`.
#'
#' @keywords internal
#'
#' @md

validate_confidence_level <- function(confidence_level = 0.95){

  #VALIDATORS:
  #Check for values >= 0, >=100 or ~= 1
  ifelse(
    #Check that confidence_level is a numeric variable
    (length(confidence_level) != 1),

    cli::cli_abort("More than one `confidence_level` provided: {confidence_level}"),

  ifelse(
    #Check that confidence_level is a numeric variable
    (!is.numeric(confidence_level) | is.na(confidence_level)),

      cli::cli_abort("Non-numeric `confidence_level` provided: {confidence_level}"),

  ifelse(
    #Check for the case of invalid confidence intervals
    (confidence_level < 0 | confidence_level > 100),

      cli::cli_abort("Confidence level of {confidence_level} is not valid.
                     Choose a confidence level between 0 and 1."),
  ifelse(
    #Check for the case of confidence intervals indistinguishable from 1
    (confidence_level <= 1 & abs(confidence_level - 1) < 1.e-6) |
      (confidence_level >= 1 & abs(confidence_level/100 - 1) < 1.e-6),

      cli::cli_abort("Confidence level of {confidence_level} is indistinguishable from 1 and
                           will not give accurate results."),

  ifelse(

    #Check for the case of confidence intervals indistinguishable from 0
    (confidence_level <= 1 & abs(confidence_level - 0) < 1.e-6) |
      (confidence_level >= 1 & abs(confidence_level/100 - 0) < 1.e-6),

      cli::cli_abort("Confidence level of {confidence_level} is indistinguishable from 0 and
                           will not give accurate results."),

    TRUE

  )))))

  #If value > 1 then divide
  if (confidence_level > 1){
    .confidence_level <- confidence_level/100
  } else {
    .confidence_level <- confidence_level
  }

  return(.confidence_level)

}

#' @title Validate number of cores
#'
#' @description
#' Internal function to validate the number of cores variable in the model
#'
#' @inheritParams pif_survey_bootstrap
#'
#' @returns An integer with the number of cores to use.
#' @importFrom parallel detectCores
#' @keywords internal
validate_number_of_cores <- function(num_cores){

  ifelse(
    #Check that cores is not a vector
    (length(num_cores) > 1),

    cli::cli_abort("Provided more than one number of cores, `num_cores`: {num_cores}."),

  ifelse(
    #Check that cores is a numeric variable
    (!is.numeric(num_cores) & !is.null(num_cores)),

    cli::cli_abort("Non-numeric number of cores, `num_cores` provided: {num_cores}."),

  ifelse(
    #Check negative or zero cores
    (!is.null(num_cores) & as.integer(num_cores) <= 0),

    cli::cli_abort("Invalid number of cores, `num_cores` provided: {num_cores}."),

  ifelse(
    #Check more cores than available
    (!is.null(num_cores) & num_cores > parallel::detectCores()),

    cli::cli_abort("Cores requested: {num_cores} > available cores: {parallel::detectCores()}."),

    TRUE
  ))))

  #Return either NULL or number of cores
  #https://stackoverflow.com/a/49388710/5067372
  .num_cores <- if (is.null(num_cores)) NULL else as.integer(num_cores)
  return(.num_cores)

}

#' @title Validate n_bootstrap_samples
#'
#' @description
#' Internal function to validate the number of bootstrap samples
#'
#' @param n_bootstrap_samples
#' @returns A valid `n_bootstrap_samples`
#' @keywords internal

validate_n_boostrap_samples <- function(n_bootstrap_samples){
  return(n_bootstrap_samples)
}

#' @title Validate design
#'
#' @description
#' Internal function to validate the design whether its a `svyrep.design` or `survey.design`.
#'
#' @param design
#' @returns A `svyrep.design`. If `design` was already a `svyrep.design` returns `design`. Else
#' it returns a bootstrap `svyrep.design` with as many as `n_bootstrap_samples`
#' @importFrom svrep as_bootstrap_design
#' @keywords internal
validate_survey_design <- function(design, n_bootstrap_samples, ...){

  .n_bootstrap_samples <- validate_n_boostrap_samples(n_bootstrap_samples = n_bootstrap_samples)

  #Throw warning if object is already a design and number of samples is set
  if (inherits(design, "svyrep.design")){
    if (!is.null(.n_bootstrap_samples)){
      cli::cli_alert_warning("Ignoring `n_bootstrap_samples` as design object is already a
        `svyrep.design` with {dim(weights(design, type = 'analysis'))[2]} replicates")
    }
    .design <- design

  #Create a survey replicate design if not already set
  } else if (inherits(design, "survey.design") | inherits(design, "survey.design2")) {

    .design <- svrep::as_bootstrap_design(design, replicates = .n_bootstrap_samples, ...)

  #Fail if object had no survey attributes
  } else {
    cli::cli_abort("design is not a `survey.design`, `survey.design2` or `svyrep.design`.
                   Use `survey::survey` to generate your design.")
  }

  return(.design)

}



#VALIDATE THETA DISTRIBUTION HAS A PARAMETER n
#VALIDATE NUMBER OF CORES
#VALIDATE additional_theta_arguments WHEN IS A LIST TO NOT HAVE AN N
#THE DERIVATIVE FUNCTION MUST BE A NAMED VECTOR WITH COLUMN NAMES AS
x <- c()
