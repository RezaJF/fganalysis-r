# Test cases for parallel_compute_purchase_frequencies_for_VNRs and compute_purchase_frequency
library(testthat)
library(dplyr)

# Dummy data for testing with Date columns
set.seed(123)
base_date <- as.Date("2020-01-01")
test_data <- data.frame(
  VNR = rep(c("A", "B"), each = 5),
  FINNGENID = rep(c("ID1", "ID2"), times = 5),
  APPROX_EVENT_DAY = c(
    base_date,
    base_date + 5,
    base_date + 10,
    base_date + 25,   # within gap
    base_date + 40,   # outside gap
    base_date + 1,
    base_date + 6,
    base_date + 11,
    base_date + 26,   # within gap
    base_date + 41    # outside gap
  ),
  PackageSize = rep(10, 10),
  DDDPerPack = rep(10, 10),
  Substance = rep("TestDrug", 10),
  ATC = rep("C10AA01", 10)
)

# Test compute_purchase_frequency

test_that("compute_purchase_frequency returns correct intervals with Date column", {
  res <- compute_purchase_frequency(test_data %>% filter(VNR == "A"), gap = 15)
  expect_true(is.data.frame(res))
  expect_true(all(res$VNR == "A"))
  expect_true(all(res$cadence >= 0))
  # Only adjacent purchases within PackageSize + gap (10 + 15 = 25 days) are returned
  res_id2 <- compute_purchase_frequency(test_data %>% filter(VNR == "A", FINNGENID == "ID2"), gap = 15)
  expect_equal(res_id2$cadence, c(20))
  # Purchases at day 0, 5, 10, 25, 40 (for ID1) should yield intervals of 5, 5
  res_id1 <- compute_purchase_frequency(test_data %>% filter(VNR == "A", FINNGENID == "ID1"), gap = 15)
  expect_equal(res_id1$cadence, c(10))
})

# Test parallel_compute_purchase_frequencies_for_VNRs

test_that("parallel_compute_purchase_frequencies_for_VNRs returns intervals for all VNRs with Date column", {
  res <- parallel_compute_purchase_frequencies_for_VNRs(test_data, gap = 15, n_workers = 1)
  expect_true(is.data.frame(res))
  expect_true(all(res$VNR %in% c("A", "B")))
  expect_true(all(res$cadence >= 0))
  # Only adjacent purchases within PackageSize + gap (10 + 15 = 25 days) are returned for both VNRs
  res_a <- res %>% filter(VNR == "A", FINNGENID == "ID2")
  expect_equal(res_a$cadence, c(20))
  res_b <- res %>% filter(VNR == "B", FINNGENID == "ID2")
  expect_equal(res_b$cadence, c(10))
})

test_that("parallel_compute_purchase_frequencies_for_VNRs works with multiple workers and Date column", {
  res <- parallel_compute_purchase_frequencies_for_VNRs(test_data, gap = 15, n_workers = 2)
  expect_true(is.data.frame(res))
  expect_true(all(res$VNR %in% c("A", "B")))
  res_a <- res %>% filter(VNR == "A", FINNGENID == "ID2")
  expect_equal(res_a$cadence, c(20))
})
