.onLoad <- function(...) {
  S7::methods_register()
  cli::cli_h1("pifpaf")
  cli::cli_text("{gsub('.\n\nA BibTeX entry.*', '', format(citation('pifpaf'))[2])}")
  cli::cli_rule()

  if (getRversion() < "4.3.0"){
    cli::cli_alert_danger("Please update your `R` version to a version > 4.3.0")
  }
}
