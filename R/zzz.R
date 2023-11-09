.onLoad <- function(...) {
  S7::methods_register()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0("pifpaf ",
           paste0(rep("\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501", 12), collapse = ""),
           "\u2501",
           paste0("\n", gsub('.\n\nA BibTeX entry.*', '',
                             base::format(utils::citation('pifpaf'))[2]),"\n"),
           paste0(rep("\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501", 13), collapse = ""))
  )
}

