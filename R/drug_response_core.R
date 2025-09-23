#' @import dplyr
#' @importFrom stats median
#' @importFrom dplyr bind_rows select mutate filter left_join group_by summarise
#' @importFrom utils head
#' @importFrom rlang .data
NULL

# Declare global variables to avoid R CMD check notes
utils::globalVariables(c("FINNGENID", "EVENT_AGE", "first_drug_age", "time_to_drug"))

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
  class(obj) <- "drug.response"
  obj
}

#' @title Create drug response object
#' @param conn A `fg_data_connection` object containing the data sources
#' @param lablist A character vector of OMOP concept IDs for the labs of interest
#' @param druglist A character vector of ATC drug codes
#' @param before_period A numeric vector of length 2 defining the before period (e.g., c(-1, 0) for 1 year to 0 before drug)
#' @param after_period A numeric vector of length 2 defining the after period (e.g., c(0.1, 1) for 0.1 to 1 year after drug)
#' @param finngen_ids Optional character vector of FINNGENIDs to restrict analysis
#' @param remove_outliers_sd Optional numeric value (1-6) to remove outliers based on standard deviations
#' @return drug.response object
#' @export
create_drug_response <- function(conn, lablist, druglist,
                                 before_period, after_period,
                                 finngen_ids = NULL, remove_outliers_sd = NULL) {

  # Validate all input parameters immediately
  if (!inherits(conn, "fg_data_connection")) {
    stop("conn must be an fg_data_connection object")
  }

  if (is.null(lablist) || is.null(druglist) || is.null(before_period) || is.null(after_period)) {
    stop("lablist, druglist, before_period, and after_period are required parameters")
  }

  # Validate parameter types and formats
  if (!is.character(lablist) || length(lablist) == 0) {
    stop("lablist must be a non-empty character vector")
  }

  if (!is.character(druglist) || length(druglist) == 0) {
    stop("druglist must be a non-empty character vector")
  }

  if (!is.numeric(before_period) || length(before_period) != 2) {
    stop("before_period must be a numeric vector of length 2")
  }

  if (!is.numeric(after_period) || length(after_period) != 2) {
    stop("after_period must be a numeric vector of length 2")
  }

  # Validate optional parameters
  if (!is.null(remove_outliers_sd)) {
    if (!is.numeric(remove_outliers_sd) || remove_outliers_sd < 1 || remove_outliers_sd > 6) {
      stop("remove_outliers_sd must be a numeric value between 1 and 6")
    }
  }

  if (!is.null(finngen_ids) && (!is.character(finngen_ids) || length(finngen_ids) == 0)) {
    stop("finngen_ids must be a non-empty character vector or NULL")
  }

  # Extract data from connection object
  kanta <- conn$labs
  phenos <- conn$pheno

  print("Querying lab measurements...")
  lab_measurements <- get_lab_measurements(all_labs = kanta,
                                           lablist = lablist,
                                           require_values = TRUE,
                                           finngen_ids = finngen_ids)

  if (!is.null(remove_outliers_sd)) {

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

  # get_drug_purchases returns ATC column by default (renamed from CODE1)
  dr_first_purchase <- drug_purchases %>% group_by(.data$FINNGENID) %>% arrange(.data$EVENT_AGE) %>%
    dplyr::summarize(n = n(), first_drug_age = first(.data$EVENT_AGE), first_drug=first(.data$ATC))

  lab_measurements <- dplyr::left_join(lab_measurements, dr_first_purchase, by = "FINNGENID")
  lab_measurements <- lab_measurements %>% mutate(time_to_drug = .data$first_drug_age - .data$EVENT_AGE)

  print("generating response summary...")

  lab_response <- generate_response_summary(lab_measurements, before_period, after_period)
  cat("Number of individuals with response data: ", nrow(lab_response), "\n")

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

  # Determine lab period relative to the drug initiation
  # IMPORTANT: time_to_drug = first_drug_age - EVENT_AGE
  # So: positive time_to_drug = BEFORE drug; negative time_to_drug = AFTER drug
  # But the period parameters use the opposite convention for compatibility:
  # before_period uses negative values (e.g., c(-1, 0) for 1 year to 0 before drug)
  # after_period uses positive values (e.g., c(0.1, 1) for 0.1 to 1 year after drug)
  # We need to flip the signs when matching
  lab_measurements <- lab_measurements %>% mutate(lab_period = case_when(
    dplyr::between(.data$time_to_drug, -before_period[2], -before_period[1]) ~ 'Before',
    dplyr::between(.data$time_to_drug, -after_period[2], -after_period[1]) ~ 'After',
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

#' @title Get Lab Measurements Before First Drug Purchase
#' @description Retrieves lab measurements for a given set of labs, filtering for a specific time window
#' before the first purchase of a specified drug. This function includes all measurements for individuals
#' who have never purchased the drug. It is a standalone function that handles data retrieval, filtering,
#' and outlier removal. The output includes a column `n_measurements` showing the total number of
#' measurements per individual.
#' @param conn A `fg_data_connection` object.
#' @param lablist A character vector of OMOP concept IDs for the labs of interest.
#' @param druglist A character vector of ATC drug codes.
#' @param months_before A numeric value specifying the time window in months before the first drug
#' purchase to include lab measurements. Defaults to 3.
#' @param remove_outliers_sd An optional numeric value specifying the number of standard deviations
#' to use for outlier removal. Values outside `mean ± sd * remove_outliers_sd` will be removed.
#' @param winsorize_pct An optional numeric value between 0 and 0.5 for Winsorizing the lab values.
#' Represents the percentage to winsorize on each tail. For example, 0.05 (5%) will cap values
#' below the 5th percentile and above the 95th percentile.
#' @param range_sd_filter An optional list with `lower_bound`, `upper_bound`, and `nsd` for
#' range-based standard deviation filtering.
#' @return A data frame of lab measurements with an `n_measurements` column, compatible with `calculate_blup_slopes`.
#' @export
get_measurements_before_drug <- function(conn, lablist, druglist, months_before,
                                         remove_outliers_sd = NULL, winsorize_pct = NULL,
                                         range_sd_filter = NULL) {

  if (!is.null(remove_outliers_sd) && !is.null(winsorize_pct)) {
    stop("Please specify only one outlier removal method: `remove_outliers_sd` or `winsorize_pct`.")
  }
  if ((!is.null(remove_outliers_sd) || !is.null(winsorize_pct)) && !is.null(range_sd_filter)) {
    stop("`range_sd_filter` cannot be used with `remove_outliers_sd` or `winsorize_pct`.")
  }

  # 1. Get all relevant lab measurements and drug purchases
  lab_measurements <- get_lab_measurements(conn$labs, lablist, require_values = TRUE)

  first_purchases <- get_first_purchase(conn$pheno, druglist) %>%
    select(FINNGENID, first_drug_age = EVENT_AGE)

  # 2. Join lab data with first purchase data
  all_measurements <- left_join(lab_measurements, first_purchases, by = "FINNGENID") %>%
    mutate(time_to_drug = .data$first_drug_age - .data$EVENT_AGE)

  # 3. Handle outlier removal
  if (!is.null(range_sd_filter)) {
    # Validate the filter parameters
    required_params <- c("lower_bound", "upper_bound", "nsd")
    if (!is.list(range_sd_filter) || !all(required_params %in% names(range_sd_filter))) {
      stop("`range_sd_filter` must be a list containing `lower_bound`, `upper_bound`, and `nsd`.")
    }

    # Calculate mean and sd on values within the specified range
    ranged_data <- all_measurements %>%
      filter(.data$MEASUREMENT_VALUE_HARMONIZED >= range_sd_filter$lower_bound &
             .data$MEASUREMENT_VALUE_HARMONIZED <= range_sd_filter$upper_bound)

    mean_val <- mean(ranged_data$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
    sd_val <- sd(ranged_data$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)

    # Define outlier bounds based on the ranged statistics
    lower_bound_stat <- mean_val - (range_sd_filter$nsd * sd_val)
    upper_bound_stat <- mean_val + (range_sd_filter$nsd * sd_val)

    original_rows <- nrow(all_measurements)
    all_measurements <- all_measurements %>%
      filter(.data$MEASUREMENT_VALUE_HARMONIZED >= lower_bound_stat & .data$MEASUREMENT_VALUE_HARMONIZED <= upper_bound_stat)
    removed_rows <- original_rows - nrow(all_measurements)
    print(paste("Removed", removed_rows, "outliers using range_sd_filter."))
  }

  if (!is.null(remove_outliers_sd)) {
    mean_val <- mean(all_measurements$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
    sd_val <- sd(all_measurements$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
    lower_bound <- mean_val - (remove_outliers_sd * sd_val)
    upper_bound <- mean_val + (remove_outliers_sd * sd_val)

    all_measurements <- all_measurements %>%
      filter(.data$MEASUREMENT_VALUE_HARMONIZED >= lower_bound & .data$MEASUREMENT_VALUE_HARMONIZED <= upper_bound)
  }

  if (!is.null(winsorize_pct)) {
    all_measurements <- all_measurements %>%
      mutate(MEASUREMENT_VALUE_HARMONIZED = winsorize_vector(.data$MEASUREMENT_VALUE_HARMONIZED, winsorize_pct))
  }

  # 4. Filter measurements based on the time window for exposed individuals
  time_window_years <- months_before / 12

  exposed_measurements <- all_measurements %>%
    filter(!is.na(.data$first_drug_age) & .data$time_to_drug >= 0 & .data$time_to_drug >= time_window_years)

  unexposed_measurements <- all_measurements %>%
    filter(is.na(.data$first_drug_age))

  # 5. Combine and add measurement count per individual
  final_measurements <- bind_rows(exposed_measurements, unexposed_measurements)

  # Add total measurement count per individual
  measurement_counts <- final_measurements %>%
    group_by(.data$FINNGENID) %>%
    summarise(n_measurements = n(), .groups = "drop")

  final_measurements <- final_measurements %>%
    left_join(measurement_counts, by = "FINNGENID")

  return(final_measurements)
}

#' @title Get Median Lab Values Before Drug Purchase with MAD Outlier Removal
#' @description Calculates the median lab value for each individual within a specified time window before
#' the first purchase of a given drug. This function includes robust outlier removal using the
#' Median Absolute Deviation (MAD) method and generates diagnostic plots.
#' @param conn A `fg_data_connection` object.
#' @param lablist A character vector of OMOP concept IDs.
#' @param druglist A character vector of ATC drug codes.
#' @param months_before The time window in months before the first drug purchase (default: 1).
#' @param remove_outliers_mad_th The threshold for MAD-based outlier removal (default: 5).
#' @param generate_plots Logical, whether to generate and save diagnostic plots (default: `FALSE`).
#' @param output_dir The directory to save outputs (default: `"."`).
#' @param output_file_prefix A prefix for output file names.
#' @return A data frame with median lab values per individual, ready for GWAS analysis.
#' @export
get_median_pre_drug <- function(conn, lablist, druglist, months_before = 1,
                                remove_outliers_mad_th = 5,
                                generate_plots = FALSE,
                                output_dir = ".",
                                output_file_prefix = "") {

  # 1. Get measurements before drug purchase
  measurements <- get_measurements_before_drug(
    conn = conn,
    lablist = lablist,
    druglist = druglist,
    months_before = months_before
  )

  # 2. MAD outlier removal
  measurements_mad <- measurements
  if (!is.null(remove_outliers_mad_th)) {
    measurements_mad <- measurements %>%
      filter(.data$MEASUREMENT_VALUE_HARMONIZED %in% filter_outliers_mad(.data$MEASUREMENT_VALUE_HARMONIZED, th = remove_outliers_mad_th))
  }

  # 3. Generate diagnostic plots if requested
  if (generate_plots) {

    # Distribution plot before and after MAD removal
    data_before <- data.frame(value = measurements$MEASUREMENT_VALUE_HARMONIZED, group = "Before")
    data_after <- data.frame(value = measurements_mad$MEASUREMENT_VALUE_HARMONIZED, group = "After")
    plot_data <- rbind(data_before, data_after)

    dist_plot <- gghistogram(plot_data, x = "value",
      add = "mean", rug = TRUE,
      color = "group", fill = "group",
      palette = c("#00AFBB", "#E7B800"),
      title = "Distribution Before and After MAD Outlier Removal",
      xlab = "Lab Value",
      ylab = "Density"
    )
    ggsave(file.path(output_dir, paste0(output_file_prefix, "_mad_distribution.png")), plot = dist_plot)

    # Violin plot of median distribution between males and females
    sex_col <- if("SEX_IMPUTED" %in% covariate_cols) "SEX_IMPUTED" else if("SEX" %in% covariate_cols) "SEX" else NULL

    if (!is.null(sex_col)) {
      # Prepare sex column for visualization
      if (sex_col == "SEX_IMPUTED") {
        measurements <- measurements %>%
          mutate(SEX_VIS = case_when(
            .data[[sex_col]] == 0 ~ "Male",
            .data[[sex_col]] == 1 ~ "Female",
            TRUE ~ "Unknown"
          ))
        measurements_mad <- measurements_mad %>%
          mutate(SEX_VIS = case_when(
            .data[[sex_col]] == 0 ~ "Male",
            .data[[sex_col]] == 1 ~ "Female",
            TRUE ~ "Unknown"
          ))
      } else {
        measurements <- measurements %>%
          mutate(SEX_VIS = case_when(
            toupper(.data[[sex_col]]) %in% c("M", "MALE", "1") ~ "Male",
            toupper(.data[[sex_col]]) %in% c("F", "FEMALE", "2") ~ "Female",
            TRUE ~ "Unknown"
          ))
        measurements_mad <- measurements_mad %>%
          mutate(SEX_VIS = case_when(
            toupper(.data[[sex_col]]) %in% c("M", "MALE", "1") ~ "Male",
            toupper(.data[[sex_col]]) %in% c("F", "FEMALE", "2") ~ "Female",
            TRUE ~ "Unknown"
          ))
      }

      pre_mad_sex <- measurements %>%
        filter(.data$SEX_VIS != "Unknown") %>%
        group_by(.data$FINNGENID, .data$SEX_VIS) %>%
        summarise(median_val = median(.data$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE), .groups = "drop") %>%
        mutate(group = "Before MAD")

      post_mad_sex <- measurements_mad %>%
        filter(.data$SEX_VIS != "Unknown") %>%
        group_by(.data$FINNGENID, .data$SEX_VIS) %>%
        summarise(median_val = median(.data$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE), .groups = "drop") %>%
        mutate(group = "After MAD")

      violin_data <- rbind(pre_mad_sex, post_mad_sex)

      violin_plot <- ggviolin(violin_data, x = "group", y = "median_val", fill = "group",
               palette = c("#FC4E07", "#00AFBB"),
               add = "boxplot", add.params = list(fill = "white")) +
        facet_wrap(~SEX_VIS) +
        stat_compare_means(label = "p.signif") +
        labs(title = "Median Lab Values by Sex (Before and After MAD)",
             x = "Group", y = "Median Lab Value")

      ggsave(file.path(output_dir, paste0(output_file_prefix, "_sex_violin.png")), plot = violin_plot)
    }
  }

  # 4. Calculate median values
  median_values <- measurements_mad %>%
    group_by(FINNGENID) %>%
    summarise(
      median_value = median(.data$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE),
      .groups = "drop"
    )

  # 5. Format for output
  output_df <- data.frame(
    FID = median_values$FINNGENID,
    IID = median_values$FINNGENID,
    median_lab_value = median_values$median_value
  )
  colnames(output_df)[3] <- paste0(lablist[1], "_median")

  # 6. Save to file
  output_file <- file.path(output_dir, paste0(output_file_prefix, "_", lablist[1], "_DF13_median.tsv"))
  write.table(output_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

  return(output_df)
}