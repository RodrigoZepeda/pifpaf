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

validate_confidence_level <- function(confidence_level = 0.95) {
  # VALIDATORS:
  # Check for values >= 0, >=100 or ~= 1
  ifelse(
    # Check that confidence_level is a numeric variable
    (length(confidence_level) != 1),
    cli::cli_abort("More than one `confidence_level` provided: {confidence_level}"),
    ifelse(
      # Check that confidence_level is a numeric variable
      (!is.numeric(confidence_level) | is.na(confidence_level)),
      cli::cli_abort("Non-numeric `confidence_level` provided: {confidence_level}"),
      ifelse(
        # Check for the case of invalid confidence intervals
        (confidence_level < 0 | confidence_level > 100),
        cli::cli_abort("Confidence level of {confidence_level} is not valid.
                     Choose a confidence level between 0 and 1."),
        ifelse(
          # Check for the case of confidence intervals indistinguishable from 1
          (confidence_level <= 1 & abs(confidence_level - 1) < 1.e-6) |
            (confidence_level >= 1 & abs(confidence_level / 100 - 1) < 1.e-6),
          cli::cli_abort("Confidence level of {confidence_level} is indistinguishable from 1 and
                           will not give accurate results."),
          ifelse(

            # Check for the case of confidence intervals indistinguishable from 0
            (confidence_level <= 1 & abs(confidence_level - 0) < 1.e-6) |
              (confidence_level >= 1 & abs(confidence_level / 100 - 0) < 1.e-6),
            cli::cli_abort("Confidence level of {confidence_level} is indistinguishable from 0 and
                           will not give accurate results."),
            TRUE
          )
        )
      )
    )
  )

  # If value > 1 then divide
  if (confidence_level > 1) {
    confidence_level <- confidence_level / 100
  } else {
    confidence_level <- confidence_level
  }

  return(confidence_level)
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
validate_number_of_cores <- function(num_cores) {
  ifelse(
    # Check that cores is not a vector
    (length(num_cores) > 1),
    cli::cli_abort("Provided more than one number of cores, `num_cores`: {num_cores}."),
    ifelse(
      # Check that cores is a numeric variable
      (!is.numeric(num_cores) & !is.null(num_cores)),
      cli::cli_abort("Non-numeric number of cores, `num_cores` provided: {num_cores}."),
      ifelse(
        # Check negative or zero cores
        (!is.null(num_cores) & as.integer(num_cores) <= 0),
        cli::cli_abort("Invalid number of cores, `num_cores` provided: {num_cores}."),
        ifelse(
          # Check more cores than available
          (!is.null(num_cores) & num_cores > parallel::detectCores()),
          cli::cli_abort("Cores requested: {num_cores} > available cores: {parallel::detectCores()}."),
          TRUE
        )
      )
    )
  )

  # Return either NULL or number of cores
  # https://stackoverflow.com/a/49388710/5067372
  num_cores <- if (is.null(num_cores)) NULL else as.integer(num_cores)
  return(num_cores)
}

#' @title Validate n_bootstrap_samples
#'
#' @description
#' Internal function to validate the number of bootstrap samples
#'
#' @param n_bootstrap_samples The number of samples for the boostrap
#' @returns A valid `n_bootstrap_samples` (integer)
#' @keywords internal
validate_n_boostrap_samples <- function(n_bootstrap_samples) {
  ifelse(
    (!is.null(n_bootstrap_samples) && !is.numeric(n_bootstrap_samples)),
    cli::cli_abort("Non-numeric variable in `n_bootstrap_samples`: {n_bootstrap_samples}"),
    ifelse(
      (length(n_bootstrap_samples) > 1),
      cli::cli_abort("More than one number of bootstrap samples provided: {n_bootstrap_samples}"),
      ifelse(
        (!is.null(n_bootstrap_samples) && as.integer(n_bootstrap_samples) <= 0),
        cli::cli_abort("Invalid number of samples for bootstrap:
                     `n_bootstrap_samples` = {n_bootstrap_samples} <= 0"),
        TRUE
      )
    )
  )

  n_bootstrap_samples <- if (is.null(n_bootstrap_samples)) NULL else as.integer(n_bootstrap_samples)
  return(n_bootstrap_samples)
}

#' @title Validate design
#'
#' @description
#' Internal function to validate the design whether its a `svyrep.design` or `survey.design`.
#'
#' @param design A `svyrep.design`, `survey.design` or `survey.design2` object.
#'
#' @returns A `svyrep.design`. If `design` was already a `svyrep.design` returns `design`. Else
#' it returns a bootstrap `svyrep.design` with as many as `n_bootstrap_samples`
#' @importFrom svrep as_bootstrap_design
#' @keywords internal
validate_survey_design <- function(design, n_bootstrap_samples, ...) {
  n_bootstrap_samples <- validate_n_boostrap_samples(n_bootstrap_samples = n_bootstrap_samples)

  # Throw warning if object is already a design and number of samples is set
  if (inherits(design, "svyrep.design")) {
    if (!is.null(n_bootstrap_samples)) {
      cli::cli_alert_warning("Ignoring `n_bootstrap_samples` as design object is already a
        `svyrep.design` with {dim(weights(design, type = 'analysis'))[2]} replicates")
    }
    design <- design

    # Create a survey replicate design if not already set
  } else if (inherits(design, "survey.design") | inherits(design, "survey.design2")) {
    if (is.null(n_bootstrap_samples)) {
      n_bootstrap_samples <- 2000
      cli::cli_alert_warning("Setting default number of bootstrap samples to:
                             {`n_bootstrap_samples`}. Consider a larger number
                             for reporting.")
    }
    design <- svrep::as_bootstrap_design(design, replicates = n_bootstrap_samples, ...)

    # Fail if object had no survey attributes
  } else {
    cli::cli_abort("design is not a `survey.design`, `survey.design2` or `svyrep.design`.
                   Use `survey::svydesign` to generate your design.")
  }

  return(design)
}

#' @title Validate parallel setup
#'
#' @description
#' Internal function to validate the parallelization and return the adecuate %do% operator
#' for foreach.
#'
#' @param parallel Boolean indicating whether to run argument in parallel
#' @param num_cores Number of cores to run the parallelization for
#' @returns A function to parallelize [foreach::foreach()] either `%do%` or `%dopar%`.
#' @keywords internal

validate_parallel_setup <- function(parallel, num_cores) {
  do_command <- ifelse(
    (!is.logical(parallel)),
    cli::cli_abort("Non-logical argument to `parallel` set it to `TRUE` or `FALSE`"),
    ifelse(
      (parallel & num_cores > 1),
      foreach::`%dopar%`,
      foreach::`%do%`
    )
  )

  return(do_command)
}

#' @title Validate the theta_distribution function
#'
#' @description
#' A function that validates that the `theta_distribution` is indeed a function
#' and that one of its arguments is `n`.
#'
#' @param `theta_distribution` the random number simulator for theta
#' @param `additional_theta_arguments` list of additional arguments to use in `theta_distribution`
#'
#' @return A random number generating function to simulate theta from.
#' @keywords internal
validate_theta_arguments <- function(theta_distribution, additional_theta_arguments, theta) {
  if (!is.function(theta_distribution) && theta_distribution == "default") {
    ifelse(
      (length(additional_theta_arguments) == 1),

      # Sigma as a variable
      theta_args <- list(
        "sigma" = as.matrix(additional_theta_arguments[[1]]),
        "mean" = theta
      ),
      ifelse(
        is.vector(additional_theta_arguments),

        # Diagonal matrix for independent
        theta_args <- list(
          "sigma" = diag(additional_theta_arguments),
          "mean" = theta
        ),

        # Default (matrix)
        theta_args <- list(
          "sigma" = additional_theta_arguments,
          "mean" = theta
        )
      )
    )
  } else {
    ifelse(
      ("n" %in% names(additional_theta_arguments)),
      cli::cli_abort("Cannot have an argument named `n` in `additional_theta_arguments` as it
                     collides with the number of simulations argument automatically set up
                     by this program."),
      theta_args <- additional_theta_arguments
    )
  }

  return(theta_args)
}

#' @title Validate the theta_distribution function
#'
#' @description
#' A function that validates that the `theta_distribution` is indeed a function
#' and that one of its arguments is `n`.
#'
#' @param `theta_distribution` the random number simulator for theta
#' @param `additional_theta_arguments` list of additional arguments to use in `theta_distribution`
#'
#' @return A random number generating function to simulate theta from.
#' @keywords internal
validate_theta_distribution <- function(theta_distribution, additional_theta_arguments) {
  ifelse(
    (!is.function(theta_distribution) && theta_distribution != "default"),
    cli::cli_abort("`theta_distribution` is not a function nor value 'default'. Choose
                   'default'  for multivariate Gaussian or use a function that generates
                   random samples for theta"),
    ifelse(
      (is.function(theta_distribution) && !("n" %in% methods::formalArgs(theta_distribution))),
      cli::cli_abort("`theta_distribution` must have a parameter `n` for generating a random sample
                   of size `n`. If you are using a function that doesn't contain that try
                   reparametrizing the function:
                   `theta_distribution = function(n, other arguments){{
                        original_function(nsim = n, other arguments)
                   }}`"),
      ifelse(
        (is.function(theta_distribution) &&
          !all(names(additional_theta_arguments) %in% methods::formalArgs(theta_distribution))),
        {
          unknown_args <- !(names(additional_theta_arguments) %in%
            methods::formalArgs(theta_distribution))

          cli::cli_alert_warning("Arguments `{names(additional_theta_arguments)[unknown_args]}` are
                             not present in `theta_distribution`")

          theta_distribution <- theta_distribution
        },
        ifelse(
          (!is.function(theta_distribution) && theta_distribution == "default"),
          theta_distribution <- mvtnorm::rmvnorm,
          TRUE
        )
      )
    )
  )

  return(theta_distribution)
}

#' @title Validate the is_paf parameter
#'
#' @description
#' A function that validates that the `is_paf` argument is valid
#'
#' @param `is_paf` boolean determining if `pif` estimation is actually `paf`
#'
#' @return A boolean indicating whether is paf or pif
#' @keywords internal
validate_is_paf <- function(is_paf) {
  if (!is.logical(is_paf)) {
    cli::cli_abort("Invalid `is_paf`: {is_paf}. Set to `FALSE` or `TRUE` if it corresponds
                   to a Population Attributable Fraction.")
  }
  return(is_paf)
}


#' @title Validate the relative risk function
#'
#' @description
#' A function that validates that the `rr` argument is valid
#'
#' @param `rr` function
#'
#' @return The relative risk function `rr`
#' @keywords internal
validate_rr <- function(rr) {
  if ((!is.function(rr) && !is.list(rr)) || (is.list(rr) && !is.function(rr[[1]]))) {
    cli::cli_abort("Invalid `rr`: {rr}. Use a function for the relative risk.")
  }

  if (!is.list(rr)){
    rr <- list(rr)
  }

  for (i in 1:length(rr)){
    if (!all(c("X", "theta") %in% methods::formalArgs(rr[[i]]))){
      arg_not <- c("X", "theta")[!(c("X", "theta") %in% methods::formalArgs(rr[[i]]))]
      cli::cli_abort("The relative risk function `rr` should have the following arguments:
                     `X` (capitalized) and `theta`: `rr(X, theta)`.")
    }
  }

  if (is.null(names(rr))){
    names(rr) <- paste0("Relative_Risk_", 1:length(rr))
  }

  return(rr)
}

#' @title Validate the counterfactual function
#'
#' @description
#' A function that validates that the `cft` argument is valid
#'
#' @param `rr` function
#'
#' @return The counterfactual function `cft`
#' @keywords internal
validate_cft <- function(cft, is_paf) {

  if (!is.function(cft) && !is.null(cft) && (is.list(cft) && !is.function(cft[[1]]))) {
    cli::cli_abort("Invalid `cft`: {cft}. Use a function for the counterfactual.")
  }

  if (!is_paf && is.null(cft)){
    cli::cli_abort("Set up a counterfactual function `cft` or change the
                   indicator `is_paf` to `TRUE` if calculating population
                   attributable fractions `paf`")
  }

  if (!is.list(cft)){
    cft <- list(cft)
  }

  cft_names <- names(cft)
  for (i in 1:length(cft)){
    if (!is.null(cft[[i]]) && !(c("X") %in% methods::formalArgs(cft[[i]]))){
      cli::cli_abort("The counterfactual function `cft` should have an argument named `X`
                     for the data.frame of exposure.")
    }

    if (is.null(cft_names[i]) && !is.null(cft[[i]])){
      names(cft)[i] <- paste0("Counterfactual_", i)
    } else if (is.null(cft_names[i]) && is.null(cft[[i]])){
      names(cft)[i] <- "Theoretical_minimum_risk_level"
    }
  }

  return(cft)
}

#' @title Validate list of potential impact fractions
#'
#' @description
#' Internal function to validate a list of potential impact fractions
#'
#' @inheritParams pif_combine
#'
#' @returns Boolean variable indicating whether input is valid
#' @keywords internal
validate_pifs <- function(...){
  pifs     <- list(...)

  #Check they are the correct class
  validate <- sapply(pifs, function(x) inherits(x, "pif_class"))
  validate <- all(validate)
  if (validate){
    #Check the number of simulations is the same for all
    validate <- length(unique(sapply(pifs, n_bootstrap))) == 1 & validate
  }
  return(validate)
}
