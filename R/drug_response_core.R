#' @import dplyr
#' @importFrom stats median
NULL

#' @title data object returned from drug response analyse (create_drug_response)
#' @param responses data frame with response data
#' @param lab_measurements data frame with all lab measurements
#' @param drug_purchases data frame with all drug purchases
#' @param before_period vector with two elements, start and end of the before period # nolint: line_length_linter.
#' @param after_period vector with two elements, start and end of the after period # nolint: line_length_linter.
#' @return object of class drug.response
#' @export

drug.response <- function(responses, lab_measurements,
                          drug_purchases,
                          before_period, after_period) {
  obj <- list(responses = responses,
              all_measurements = lab_measurements,
              all_drug_purchases = drug_purchases,
              lab_response_period = list(before_period = before_period,
                                         after_period = after_period))
  class(obj) <- "drug.reponse"
  obj
}

#' @title Create drug response object
#' @param kanta data frame with lab measurements
#' @param phenos data frame with drug purchases
#' @param lablist vector of lab measurement concept IDs
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*).
#' @param before_period vector with two elements, start and end of the before period
#' @param after_period vector with two elements, start and end of the after period
#' @param finngen_ids vector of FINNGENIDs to filter the data
#' @param remove_outliers_sd integer, if defined, remove outliers from the lab measurements. The value is the number of standard deviations from the mean to define an outlier (i.e. 1, 2, 3, 4, 5, 6)
#' @param covariates Optional data frame with covariate data (e.g., sex, age at death)
#' @param covariate_cols Optional character vector of column names to add from the covariates data frame
#' @return drug.response object
#' @export
create_drug_response <- function(kanta, phenos, lablist, druglist,
                                 before_period, after_period, finngen_ids = NULL, remove_outliers_sd = NULL,
                                 covariates = NULL,
                                 covariate_cols = NULL) {
  print("Querying lab measurements...")
  lab_measurements <- get_lab_measurements(all_labs=kanta, lablist=lablist, finngen_ids=finngen_ids, require_values = TRUE)

  if (!is.null(remove_outliers_sd)) {
    if (!is.numeric(remove_outliers_sd) || remove_outliers_sd < 1 || remove_outliers_sd > 6) {
      stop("remove_outliers_sd must be an integer between 1 and 6")
    }

    mean_val <- mean(lab_measurements$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
    sd_val <- sd(lab_measurements$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)

    lower_bound <- mean_val - (remove_outliers_sd * sd_val)
    upper_bound <- mean_val + (remove_outliers_sd * sd_val)

    original_rows <- nrow(lab_measurements)
    lab_measurements <- lab_measurements %>% # nolint
      filter(.data$MEASUREMENT_VALUE_HARMONIZED >= lower_bound & .data$MEASUREMENT_VALUE_HARMONIZED <= upper_bound)

    removed_rows <- original_rows - nrow(lab_measurements)
    print(paste("Removed", removed_rows, "outliers based on", remove_outliers_sd, "standard deviations from the mean."))
  }

  all_fg_ids <- unique(c(lab_measurements$FINNGENID, finngen_ids))
  print("Querying purchases...")
  drug_purchases <- get_drug_purchases(phenos, druglist, all_fg_ids)

  # Get the first drug code column (either ATC if renamed, or CODE1 if not)
  drug_col <- if ("ATC" %in% colnames(drug_purchases)) "ATC" else "CODE1"

  dr_first_purchase <- drug_purchases %>% group_by(.data$FINNGENID) %>% arrange(.data$EVENT_AGE) %>%
    dplyr::summarize(n = n(), first_drug_age = first(.data$EVENT_AGE), first_drug=first(.data[[drug_col]]))

  lab_measurements <- dplyr::left_join(lab_measurements, dr_first_purchase, by = "FINNGENID")
  lab_measurements <- lab_measurements %>% mutate(time_to_drug = .data$first_drug_age - .data$EVENT_AGE)

  print("generating response summary...")

  lab_response <- generate_response_summary(lab_measurements, before_period, after_period)
  cat("Number of individuals with response data: ", nrow(lab_response), "\n")

  # Add covariates if provided
  if (!is.null(covariates) && !is.null(covariate_cols)) {
    print(paste("Adding covariates:", paste(covariate_cols, collapse = ", ")))

    # Ensure FINNGENID is in the columns for the join
    cols_to_select <- unique(c("FINNGENID", covariate_cols))

    # Check that all requested columns exist in the covariates dataframe
    missing_cols <- setdiff(cols_to_select, colnames(covariates))
    if (length(missing_cols) > 0) {
      stop(paste("The following `covariate_cols` are not in the `covariates` dataframe:",
                 paste(missing_cols, collapse = ", ")))
    }

    # Select the requested columns and ensure one row per FINNGENID
    # Collect the data if it's a lazy table
    cov_data_to_join <- covariates %>%
      select(all_of(cols_to_select)) %>%
      distinct(.data$FINNGENID, .keep_all = TRUE)

    # If covariates is a lazy table (database connection), collect it
    if (inherits(covariates, "tbl_lazy") || inherits(covariates, "tbl_sql")) {
      cov_data_to_join <- dplyr::collect(cov_data_to_join)
    }

    # Join with the final response summary data frame
    lab_response <- lab_response %>%
      left_join(cov_data_to_join, by = "FINNGENID")

    # Also join with the full measurements data frame
    lab_measurements <- lab_measurements %>%
      left_join(cov_data_to_join, by = "FINNGENID")
  }

  drug.response(responses = lab_response, lab_measurements = lab_measurements,
                drug_purchases = drug_purchases, before_period = before_period, after_period = after_period)
}


#' @title Generate response summary
#' @param lab_measurements data frame with lab measurements
#' @param before_period vector with two elements, start and end of the before period
#' @param after_period vector with two elements, start and end of the after period
#' @param summary_function function to summarize the lab measurements (default is median)
#' @return data frame with response summary
#' @export
generate_response_summary <- function(lab_measurements, before_period, after_period, summary_function=median) {

  lab_measurements <- lab_measurements %>% mutate(lab_period = case_when(
    dplyr::between(.data$time_to_drug, after_period[1], after_period[2]) ~ 'Before',
    dplyr::between(.data$time_to_drug, before_period[1], before_period[2]) ~ 'After',
    TRUE ~ NA_character_
  ))

  lab_response <- lab_measurements %>% dplyr::filter(!is.na(.data$lab_period) & !is.na(.data$MEASUREMENT_VALUE_HARMONIZED)) %>%
    dplyr::group_by(.data$FINNGENID) %>%
    dplyr::summarize(before = summary_function(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='Before'], na.rm=TRUE),
              after = summary_function(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='After'], na.rm=TRUE),
              n_before = length(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='Before']),
              n_after = length(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='After']),
              baseline_age = first(.data$first_drug_age),
              first_drug = first(.data$first_drug),
              response = .data$after - .data$before) %>% dplyr::filter(!is.na(.data$response))

  lab_response
}