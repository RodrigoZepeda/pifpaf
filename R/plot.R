#' @title Plot the population attributable fraction or the potential impact fraction
#'
#' @description
#' Creates a [ggplot2::ggplot()] object containing information on the population
#' attributable fraction or the potential impact fraction.
#'
#'
#' @param x A `pif_class` object
#' @param title The title of the plot
#' @param xaxis Element for the x axis: `counterfactual` or `relative_risk`
#' @param relative_risks Vector of relative risk names to plot (set to `NULL` to plot all)
#' @param counterfactuals Vector of counterfactual names to plot (set to `NULL` to plot all)
#' @param ... Additional (currently unused) parameters
#' @examples
#' # EXAMPLE 1
#' # Setup the survey design
#' data(ensanut)
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#'
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
#' my_pif <- pif(design,
#'   theta = log(c(1.05, 1.38, 1.21)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03, 0.025),
#'   n_bootstrap_samples = 10, parallel = FALSE
#' )
#' plot(my_pif)
#' plot(my_pif, xaxis = "counterfactual")
#' @export
#' @return A `ggplot2` object with an image of the `pif` (or `paf`)
#' @md
plot.pif_class <- function(x, ...,
                         xaxis = c("relative_risk", "counterfactual"),
                         title = NULL,
                         relative_risks = NULL,
                         counterfactuals = NULL){

  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("scales", quietly = TRUE) &&
      requireNamespace("tidyr", quietly = TRUE)){

    facet <- ifelse(xaxis[1] == "relative_risk", "counterfactual", "relative_risk")

    pif_df <- as.data.frame.pif_class(x)

    if (!is.null(relative_risks)){
      pif_df <- dplyr::filter(pif_df, !!as.symbol("relative_risk") %in% !!relative_risks)
      if (nrow(pif_df) < 1){
        cli::cli_abort("Relative risks specified in {.code relative_risks} not found. Please
                       use at least one of the following: {unique(pif_df[,'relative_risk'])}
                       or set `relative_risks = NULL` to use all.")
      }
    }

    if (!is.null(counterfactuals)){
      pif_df <- dplyr::filter(pif_df, !!as.symbol("counterfactual") %in% !!counterfactuals)
      if (nrow(pif_df) < 1){
        cli::cli_abort("Counterfactuals specified in {.code counterfactuals} not found. Please
                       use at least one of the following: {unique(pif_df[,'counterfactual'])}
                       or set `counterfactuals = NULL` to use all.")
      }
    }

    pif_df <- tidyr::pivot_wider(pif_df,
                                 names_from = "type",
                                 values_from = c(dplyr::everything(),
                                                 -dplyr::any_of(
                                                   c("counterfactual","relative_risk","type"))))
    colnames(pif_df) <- gsub("\\d|%|\\.| ", "", colnames(pif_df))

    if (is.null(title)){
      title <- get_fraction_type(x)
    }

    ggplot2::ggplot(pif_df) +
      ggplot2::geom_errorbar(
        ggplot2::aes_string(x = xaxis[1],
                            ymin = colnames(pif_df)[
                              grep("potential.*Lower|population.*Lower",colnames(pif_df))],
                            ymax = colnames(pif_df)[
                              grep("potential.*Upper|population.*Upper",colnames(pif_df))],
                            color = facet),
      ) +
      ggplot2::geom_point(
        ggplot2::aes_string(x = xaxis[1],
                            y = colnames(pif_df)[
                              grep("potential.*point|population.*point", colnames(pif_df))],
                            color = facet)
      ) +
      ggplot2::facet_wrap(stats::as.formula(paste0("~", facet))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none"
      ) +
      ggplot2::labs(
        title = title,
        y = title,
        x = ifelse(xaxis[1] == "relative_risk", "Relative Risk", "Counterfactual"),
        caption = paste("Results from", n_bootstrap(x), "bootstrap simulations")
      ) +
      ggplot2::scale_y_continuous(labels = scales::label_percent())
  } else {
    cli::cli_alert_danger(
      "Please install the {.code ggplot2}, {.code scales}, and {.code tidyr} packages to be
       able to plot: {.code install.packages(c('ggplot2','scales','tidyr'))}"
    )
  }
}


#' @title Diagnostic bootstrap plots for `pif_class`
#'
#' @description
#' Creates a [ggplot2::ggplot()] object containing diagnostic plots for the convergence
#' of the population attributable fraction and potential impact fraction
#'
#' @inheritParams plot.pif_class
#' @param level (`double`) Level for the quantiles of the plot (default to `0.95`)
#' @examples
#' # EXAMPLE 1
#' # Setup the survey design
#' data(ensanut)
#' options(survey.lonely.psu = "adjust")
#' design <- survey::svydesign(data = ensanut, ids = ~1, weights = ~weight, strata = ~strata)
#'
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
#' my_pif <- pif(design,
#'   theta = log(c(1.05, 1.38, 1.21)), rr, cft,
#'   additional_theta_arguments = c(0.01, 0.03, 0.025),
#'   n_bootstrap_samples = 10, parallel = FALSE
#' )
#' diagnostic_plot(my_pif)
#' @export
#' @return A `ggplot2` object with caterpillar plots to assess convergence
#' @md
diagnostic_plot <- function(x, relative_risks = NULL, counterfactuals = NULL, level = 0.95){

  conf_level <- validate_confidence_level(level)
  alpha_val  <- 1 - conf_level

  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("scales", quietly = TRUE) &&
      requireNamespace("cumstats", quietly = TRUE) &&
      requireNamespace("tidyr", quietly = TRUE)){

    pif_df    <- get_bootstrap_simulations(x)
    end_coef  <- coef.pif_class(x)

    if (!is.null(relative_risks)){
      pif_df   <- dplyr::filter(pif_df, !!as.symbol("relative_risk") %in% !!relative_risks)
      end_coef <- dplyr::filter(end_coef, !!as.symbol("relative_risk") %in% !!relative_risks)
      if (nrow(pif_df) < 1){
        cli::cli_abort("Relative risks specified in {.code relative_risks} not found. Please
                       use at least one of the following: {unique(pif_df[,'relative_risk'])}
                       or set `relative_risks = NULL` to use all.")
      }
    }

    if (!is.null(counterfactuals)){
      pif_df   <- dplyr::filter(pif_df, !!as.symbol("counterfactual") %in% !!counterfactuals)
      end_coef <- dplyr::filter(end_coef, !!as.symbol("counterfactual") %in% !!counterfactuals)
      if (nrow(pif_df) < 1){
        cli::cli_abort("Counterfactuals specified in {.code counterfactuals} not found. Please
                       use at least one of the following: {unique(pif_df[,'counterfactual'])}
                       or set `counterfactuals = NULL` to use all.")
      }
    }

    pif_df <- dplyr::mutate(
      dplyr::group_by_at(pif_df, c("relative_risk","counterfactual")), "sim" = 1:dplyr::n())

    pif_df <- dplyr::mutate(pif_df,
        dplyr::across(dplyr::matches("potential|population"),
           list("mean" = cumstats::cummean,
                "var" = cumstats::cumvar,
                "qlow" = ~cumstats::cumquant(., p = alpha_val/2),
                "qup" = ~cumstats::cumquant(., p = 1 - alpha_val/2))))

    pif_df  <- dplyr::ungroup(pif_df)

    #Remove 1st simulation which has na for variance
    pif_df <- tidyr::drop_na(pif_df)

    ggplot2::ggplot(pif_df) +
      ggplot2::geom_ribbon(
        ggplot2::aes_string(x = "sim",
                            ymin = colnames(pif_df)[
                              grep("potential.*qlow|population.*qlow", colnames(pif_df))],
                            ymax = colnames(pif_df)[
                              grep("potential.*qup|population.*qup", colnames(pif_df))]),
        alpha = 0.25, fill = "#41908E") +
      ggplot2::geom_line(
        ggplot2::aes_string(x = "sim",
                            y = colnames(pif_df)[
                              grep("potential.*mean|population.*mean", colnames(pif_df))]),
        color = "#992C36") +
      ggplot2::facet_grid(stats::as.formula("relative_risk ~ counterfactual")) +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Diagnostic plot",
        subtitle = get_fraction_type(x),
        x = "Simulation",
        y = "Estimation",
        caption = "Converging estimates should show stabilization toward dotted line."
      ) +
      ggplot2::geom_hline(
        ggplot2::aes_string(
          yintercept = colnames(end_coef)[
            grep("potential.*fraction|population.*fraction", colnames(end_coef))]),
        data = end_coef, linetype = "dotted", color = "gray25")
  } else {
    cli::cli_alert_danger(
      "Please install the {.code ggplot2}, {.code scales}, {.code cumstats}, and {.code tidyr}
      packages to be able to plot:
      {.code install.packages(c('ggplot2','scales','cumstats','tidyr'))}")
  }
}








