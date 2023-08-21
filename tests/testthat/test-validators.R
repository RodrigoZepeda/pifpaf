test_that("`validate_confidence_level` works", {

  #Test for the default value
  expect_equal(validate_confidence_level(), 0.95)

  #Test for NA
  expect_error(validate_confidence_level(NA_real_))

  #Test what happens if inputing a vector
  expect_error(validate_confidence_level(confidence_level = seq(0.1, 0.5, length.out = 10)))

  #Test what happens if inputing string
  expect_error(validate_confidence_level(confidence_level = "0.95"))

  #Test for confidence level > 100
  expect_error(validate_confidence_level(confidence_level = 1000))

  #Test for confidence level of 0
  expect_error(validate_confidence_level(confidence_level = 0))

  #Test for confidence level of 1
  expect_error(validate_confidence_level(confidence_level = 1))

  #Test for confidence level < 0
  expect_error(validate_confidence_level(confidence_level = -1))

  #Test for confidence level between 0 and 1
  expect_equal(validate_confidence_level(confidence_level = 0.5), 0.5)

  #Test for confidence level between 1 and 100
  expect_equal(validate_confidence_level(confidence_level = 50), 0.5)

})

test_that("`validate_number_of_cores` works", {

  #Test for the default value
  cores <- parallel::detectCores()
  expect_equal(validate_number_of_cores(num_cores = cores), cores)

  #Test for null cores
  expect_null(validate_number_of_cores(num_cores = NULL))

  #Test what happens if inputing vector
  expect_error(validate_number_of_cores(num_cores = 1:2))

  #Test what happens if inputing string
  expect_error(validate_number_of_cores(num_cores = "12"))

  #Test for number of cores > available
  expect_error(validate_number_of_cores(num_cores = cores + 1))

  #Test for negative cores and 0 cores
  expect_error(validate_number_of_cores(num_cores = 0))
  expect_error(validate_number_of_cores(num_cores = 0.2))
  expect_error(validate_number_of_cores(num_cores = -2))

  #Test for decimal cores
  expect_equal(validate_number_of_cores(num_cores = 1.3), 1)

})


