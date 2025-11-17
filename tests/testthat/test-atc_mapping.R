# Test file for ATC code mapping functionality
library(testthat)
library(fganalysis)

test_that("ATC mapping expansion works for dulaglutide", {

  # Clear cache to ensure fresh load
  clear_atc_cache()

  # Test expansion of dulaglutide current code
  expanded <- expand_atc_codes(
    atc_codes = c("A10BJ05"),  # Current dulaglutide code
    include_hierarchical = FALSE,
    verbose = FALSE
  )

  # Should include both current and historical codes
  expect_true("A10BJ05" %in% expanded)  # Current code
  expect_true("A10BX07" %in% expanded)  # Historical code

  # Test expansion of historical code
  expanded_historical <- expand_atc_codes(
    atc_codes = c("A10BX07"),  # Historical dulaglutide code
    include_hierarchical = FALSE,
    verbose = FALSE
  )

  # Should include both codes when starting from historical
  expect_true("A10BJ05" %in% expanded_historical)  # Current code
  expect_true("A10BX07" %in% expanded_historical)  # Original input
})

test_that("ATC mapping with hierarchical expansion works", {

  # Clear cache to ensure fresh load
  clear_atc_cache()

  # Test with hierarchical expansion
  expanded <- expand_atc_codes(
    atc_codes = c("A10BJ05"),
    include_hierarchical = TRUE,
    verbose = FALSE
  )

  # Should include hierarchical codes
  expect_true("A10BJ05" %in% expanded)  # Full code
  expect_true("A10BJ" %in% expanded)     # 5-character hierarchy
  expect_true("A10B" %in% expanded)      # 4-character hierarchy
  expect_true("A10" %in% expanded)       # 3-character hierarchy
})

test_that("get_atc_relationships returns correct information", {

  # Clear cache to ensure fresh load
  clear_atc_cache()

  # Get relationships for current code
  relationships <- get_atc_relationships("A10BJ05")

  expect_equal(relationships$original, "A10BJ05")
  expect_equal(relationships$current, "A10BJ05")
  expect_equal(relationships$historical, "A10BX07")
  expect_equal(relationships$description, "dulaglutide")

  # Get relationships for historical code
  relationships_old <- get_atc_relationships("A10BX07")

  expect_equal(relationships_old$original, "A10BX07")
  expect_equal(relationships_old$current, "A10BJ05")
  expect_equal(relationships_old$related, "A10BJ05")
})

test_that("Integration with get_drug_purchases works with mock data", {

  # Skip if we're not in a test environment with mock data
  skip_if_not(exists("create_mock_connection"))

  # Clear cache to ensure fresh load
  clear_atc_cache()

  # Create mock data with both old and new dulaglutide codes
  # Note: FINNGENIDs are anonymized test IDs
  pheno_data <- data.frame(
    FINNGENID = c("TEST001", "TEST002", "TEST003", "TEST004"),
    EVENT_AGE = c(50, 55, 60, 65),
    APPROX_EVENT_DAY = c(18250, 20075, 21900, 23725),
    CODE1 = c("A10BJ05", "A10BX07", "A10BJ05", "C10AA01"),  # Mix of new, old dulaglutide and a statin
    CODE2 = c("205", "205", "205", "315"),  # Reimbursement codes
    CODE3 = c("VNR123", "VNR124", "VNR125", "VNR126"),  # VNR codes
    CODE4 = c("1", "1", "1", "1"),  # Number of packs
    SOURCE = rep("PURCH", 4),
    stringsAsFactors = FALSE
  )

  labs_data <- data.frame(
    FINNGENID = c("TEST001", "TEST002", "TEST003"),
    OMOP_CONCEPT_ID = c("3004410", "3004410", "3004410"),  # HbA1c
    EVENT_AGE = c(49, 54, 59),
    VALUE = c(8.5, 9.0, 7.5),
    stringsAsFactors = FALSE
  )

  # Create mock connection
  mock_conn <- create_mock_connection(
    pheno_data = pheno_data,
    labs_data = labs_data
  )

  # Test with ATC mapping enabled (default)
  drug_purchases_mapped <- get_drug_purchases(
    conn = mock_conn,
    druglist = c("A10BJ05"),  # Search for current dulaglutide code
    use_atc_mapping = TRUE,
    use_only_reimbursement = TRUE
  )

  # Should find both TEST001 (A10BJ05) and TEST002 (A10BX07)
  expect_equal(nrow(drug_purchases_mapped), 3)  # TEST001, TEST002, and TEST003
  expect_true("TEST001" %in% drug_purchases_mapped$FINNGENID)
  expect_true("TEST002" %in% drug_purchases_mapped$FINNGENID)
  expect_true("TEST003" %in% drug_purchases_mapped$FINNGENID)
  expect_false("TEST004" %in% drug_purchases_mapped$FINNGENID)  # Different drug

  # Test with ATC mapping disabled
  drug_purchases_no_map <- get_drug_purchases(
    conn = mock_conn,
    druglist = c("A10BJ05"),
    use_atc_mapping = FALSE,
    use_only_reimbursement = TRUE
  )

  # Should only find TEST001 and TEST003 (exact match A10BJ05)
  expect_equal(nrow(drug_purchases_no_map), 2)
  expect_true("TEST001" %in% drug_purchases_no_map$FINNGENID)
  expect_false("TEST002" %in% drug_purchases_no_map$FINNGENID)  # Old code not found
  expect_true("TEST003" %in% drug_purchases_no_map$FINNGENID)
})

test_that("create_drug_response works with ATC mapping", {

  skip_if_not(exists("create_mock_connection"))

  # Clear cache to ensure fresh load
  clear_atc_cache()

  # Create comprehensive mock data
  pheno_data <- data.frame(
    FINNGENID = c(rep("TEST001", 2), rep("TEST002", 2), rep("TEST003", 1)),
    EVENT_AGE = c(50, 52, 55, 57, 60),
    APPROX_EVENT_DAY = c(18250, 18980, 20075, 20805, 21900),
    CODE1 = c("A10BJ05", "A10BJ05", "A10BX07", "A10BX07", "A10BJ05"),
    CODE2 = c("205", "205", "205", "205", "205"),
    CODE3 = paste0("VNR", 1:5),
    CODE4 = rep("1", 5),
    SOURCE = rep("PURCH", 5),
    stringsAsFactors = FALSE
  )

  labs_data <- data.frame(
    FINNGENID = c(rep("TEST001", 4), rep("TEST002", 4), rep("TEST003", 3)),
    OMOP_CONCEPT_ID = rep("3004410", 11),  # HbA1c
    EVENT_AGE = c(
      49, 49.5, 51, 53,       # TEST001: before and after
      54, 54.5, 56, 58,       # TEST002: before and after
      59, 61, 62              # TEST003: before and after
    ),
    VALUE = c(
      8.5, 8.3, 7.8, 7.5,     # TEST001: improvement
      9.0, 8.8, 8.2, 8.0,     # TEST002: improvement
      7.5, 7.0, 6.8           # TEST003: improvement
    ),
    stringsAsFactors = FALSE
  )

  mock_conn <- create_mock_connection(
    pheno_data = pheno_data,
    labs_data = labs_data
  )

  # Create drug response with ATC mapping
  response <- create_drug_response(
    conn = mock_conn,
    lablist = c("3004410"),
    druglist = c("A10BJ05"),  # Current code only
    before_period = c(-1, 0),
    after_period = c(0.1, 2),
    use_atc_mapping = TRUE  # Will find both old and new codes
  )

  # All three patients should be included
  expect_equal(nrow(response$responses), 3)

  # Check that responses are calculated correctly
  expect_true(all(response$responses$response < 0))  # All show improvement (lower HbA1c)
})

test_that("Verbose output provides clear information", {

  # Clear cache to ensure fresh load
  clear_atc_cache()

  # Capture messages
  messages <- capture.output({
    expanded <- expand_atc_codes(
      atc_codes = c("A10BJ05", "A10BX07"),
      include_hierarchical = FALSE,
      verbose = TRUE
    )
  }, type = "message")

  # Check that messages contain expected information
  expect_true(any(grepl("ATC Code Expansion", messages)))
  expect_true(any(grepl("Expanding 2 input ATC code", messages)))
  expect_true(any(grepl("expanded to", messages)))
  expect_true(any(grepl("Expansion Complete", messages)))
})
