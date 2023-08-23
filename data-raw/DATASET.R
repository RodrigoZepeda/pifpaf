## code to prepare the datasets
ensanut <- readr::read_rds(file.path("data-raw", "ensanut_sbp_sodium.rds")) |>
  dplyr::rename(age = edad) |>
  dplyr::rename(weight = ponde_f) |>
  dplyr::rename(strata = est_var) |>
  dplyr::rename(systolic_blood_pressure = mean_sist) |>
  dplyr::rename(sex = sex_name)
usethis::use_data(ensanut, overwrite = TRUE)
