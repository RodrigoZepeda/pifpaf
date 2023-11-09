test_that("`validate_confidence_level` works", {
  # Test for the default value
  expect_equal(validate_confidence_level(), 0.95)

  # Test for NA
  expect_error(validate_confidence_level(NA_real_))

  # Test what happens if inputing a vector
  expect_error(validate_confidence_level(confidence_level = seq(0.1, 0.5, length.out = 10)))

  # Test what happens if inputing string
  expect_error(validate_confidence_level(confidence_level = "0.95"))

  # Test for confidence level > 100
  expect_error(validate_confidence_level(confidence_level = 1000))

  # Test for confidence level of 0
  expect_error(validate_confidence_level(confidence_level = 0))

  # Test for confidence level of 1
  expect_error(validate_confidence_level(confidence_level = 1))

  # Test for confidence level < 0
  expect_error(validate_confidence_level(confidence_level = -1))

  # Test for confidence level between 0 and 1
  expect_equal(validate_confidence_level(confidence_level = 0.5), 0.5)

  # Test for confidence level between 1 and 100
  expect_equal(validate_confidence_level(confidence_level = 50), 0.5)
})
